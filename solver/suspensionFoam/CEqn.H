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
    
volScalarField& C = CPtr.ref();
volVectorField& Ws = WsPtr.ref();
surfaceScalarField& phip = phipPtr.ref();

Info << " - flag variables called by reference" << endl;
    
volScalarField epsEff("epsEff", turbulence->nut()/Sc);
word scheme("div(phi,C)");
word schemer("div(phir,C)");
surfaceScalarField phiC
(
    fvc::flux
    (
        phi,
        C,
        scheme
    )
  + fvc::flux
    (
        fvc::flux(Ws),
        C,
        schemer
    )
);

MULES::explicitSolve
(
    geometricOneField(),
    C,
    phi,
    phiC,
    zeroField(),
    zeroField(),
    UniformField<scalar>(0.99*0.7),
    zeroField()
);


// Apply the diffusion term separately to allow implicit solution
// and boundedness of the explicit advection

fvScalarMatrix CEqn
(
    fvm::ddt(C) - fvc::ddt(C)
  - fvm::laplacian(epsEff, C)
  ==
    fvOptions(C)
);

fvOptions.constrain(CEqn);

CEqn.solve(mesh.solver("CDiffusion"));

phiC += CEqn.flux();

C.correctBoundaryConditions();
    
fvOptions.correct(C);

Ws = settlingModelPtr->Ufall(C, Cmax, dS, rhoS, rhoF, nuF, mag(g))*g/mag(g);

Ws.correctBoundaryConditions();

// Update flux
phip = fvc::flux(U+Ws);

Info << "sediment phase volume fraction = "
        << C.weightedAverage(mesh.V()).value()
        << "  Min(C) = " << gMin(C)
        << "  Max(C) = " << gMax(C)
        << endl;
