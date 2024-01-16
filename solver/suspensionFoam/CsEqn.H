/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

volScalarField& Cs = CsPtr.ref();  // suspended load volumic fraction
volVectorField& Ws = WsPtr.ref();  // settling velocity
surfaceScalarField& phip = phipPtr.ref();  // flux(U+Ws)

// turbulent sediment diffusivity
volScalarField epsEff("epsEff", turbulence->nut()/Sc);

// MULES FORMULATION FOR SUSPENTION TRANSPORT
if (switchMULES=="on")
{
    word scheme("div(phi,Cs)");
    word schemer("div(phir,Cs)");
    surfaceScalarField phiCs
    (
        fvc::flux
        (
            phi,
            Cs,
            scheme
        )
        + fvc::flux
        (
            fvc::flux(Ws),
            Cs,
            schemer
        )
    );

    MULES::explicitSolve
    (
        geometricOneField(),
        Cs,
        phi,
        phiCs,
        zeroField(),
        zeroField(),
        UniformField<scalar>(0.99*0.7),
        zeroField()
    );

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection

    fvScalarMatrix CsEqn
    (
        fvm::ddt(Cs) - fvc::ddt(Cs)
      - fvm::laplacian(epsEff, Cs)
        ==
        fvOptions(Cs)
    );

    fvOptions.constrain(CsEqn);

    CsEqn.solve(mesh.solver("CsDiffusion"));

    phiCs += CsEqn.flux();

    Cs.correctBoundaryConditions();

    fvOptions.correct(Cs);
}

// CLASSICAL FORMULATION FOR SUSPENTION TRANSPORT
else if (switchMULES=="off")
{
    fvScalarMatrix CsEqn
    (
        fvm::ddt(Cs)
      + fvm::div(phip, Cs)
      - fvm::laplacian(epsEff, Cs)
        ==
        fvOptions(Cs)
    );

    fvOptions.constrain(CsEqn);

    CsEqn.solve();

    fvOptions.correct(Cs);
}

Ws = settlingModelPtr->Ufall(Cs, CsMax, dS, rhoS, rhoF, nuF, mag(g))*g/mag(g);

Ws.correctBoundaryConditions();

phip = fvc::flux(U + Ws);

Info << "sediment phase volume fraction = "
        << Cs.weightedAverage(mesh.V()).value()
        << "  Min(Cs) = " << gMin(Cs)
        << "  Max(Cs) = " << gMax(Cs)
        << endl;
