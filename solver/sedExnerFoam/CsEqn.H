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

// Access variables relative to suspension
volScalarField& Cs = CsPtr.ref();  // suspended load volumic fraction
volVectorField& Ws = WsPtr.ref();  // settling velocity
surfaceScalarField& phip = phipPtr.ref();  // flux(U+Ws)
surfaceScalarField& bedPhi = bedPhiPtr.ref();  // mesh flux through faces
// reference to erosion terms
volScalarField& Erosion = erosionPtr.ref();

// turbulent sediment diffusivity
volScalarField epsEff("epsEff", turbulence->nut()/Sc);
if (CsDiffusivity == "effective")
{
    volScalarField epsEff("epsEff", turbulence->nuEff()/Sc);
}

// handle erosion and deposition if bedload and suspension are switched on
if (bed.exist())
{
    // reference to faMesh
    const faMesh& aMesh = bed.aMesh();
    // Shields number
    areaVectorField& shields = shieldsPtr.ref();
    // critical Shields stress
    areaScalarField& critShields = critShieldsPtr.ref();
    

    for (const label patchID : bed.bedPatchesID())
    {
        // erosion boundary field
        scalarField& Ep = Erosion.boundaryFieldRef()[patchID];
        // reference to bed fvPatch
        const fvPatch& bedPatch = mesh.boundary()[patchID];

        // if dynamicMesh, meshPhi is needed for erosion flux correction
        if (mesh.changing())
        {
            // get meshPhi on bed patch
            bedPhi.boundaryFieldRef()[patchID] =
                mesh.setPhi().ref().boundaryFieldRef()[patchID];
        }
        const scalarField& bedPhiPatch = bedPhi.boundaryFieldRef()[patchID];
        // reference level
        scalar zbRef = 2.5 * dS.value();
        // dimensionless sediment diameter
        scalar sRatio = rhoS.value() / rhoF.value();  // density ratio
        scalar Dstar = dS.value()
            * Foam::pow(
                (sRatio - 1) * mag(g).value() / Foam::pow(nuF.value(), 2),
                1.0/3.0
            );

        forAll(bedPatch, facei)
        {
            // shear friction dimension less parameter
            scalar Tshear =
                Foam::mag(shields[facei])/critShields[facei] - 1;
            // compute reference concentration, Van Rijn 1984
            scalar cbRef =
                Foam::min(
                    0.015 * (dS.value()/zbRef)
                  * Foam::pow(Tshear * Foam::pos(Tshear), 1.5)
                  * Foam::pow(Dstar, -0.3),
                    CsMax.value()
                );
            // label of neighbour cell
            label idNCell = bedPatch.faceCells()[facei];
            // coordinates of cell centers located above bed
            const vector xCell = mesh.C()[idNCell];
            // face coordinates
            const vector xFace = aMesh.areaCentres()[facei];
            const vector nFace = aMesh.faceAreaNormals()[facei];
            const scalar Sface = aMesh.S()[facei];
            const scalar Vcell = mesh.V()[idNCell];
            scalar zi = -(xCell - xFace) & nFace;
            //scalar ci = Cs[idNCell];
            vector wsi = Ws[idNCell];
            scalar bedPhiFace = bedPhiPatch[facei];
            scalar Roi = Foam::mag(wsi) / Foam::max(epsEff[idNCell], 1e-6);
            //scalar cb = ci + cbRef * (
            //    1 - Foam::exp(-Roi * (Foam::mag(zi-zbRef))));
            // implicit, zeroGradient bc for Ws and Cs on bed needed
            // extrapolate cbRef at level
            scalar cbRefZ1 = Foam::min
                (
                    cbRef * Foam::exp(-Roi * (zi-zbRef)),
                    CsMax.value()
                );
            Ep[facei] = (cbRefZ1 / Vcell)
                * ((wsi & nFace) * Sface - bedPhiFace);
            Erosion[idNCell] = Ep[facei];
            /*Info << "dz = " << zi - zbRef
                << ", nuF = " << nuF.value()
                << ", s = " << sRatio
                << ", g = " << mag(g).value()
                << ", D* = " << Dstar
                << ", T = " << Tshear
                << ", ws = " << wsi
                << ", epsEff = " << epsEff[idNCell]
                << ", Ro = " << Roi
                << ", c1 = " << ci
                << ", cbRef = " << cbRef
                << ", cb = " << cb
                << ", cbRefZ1 = " << cbRefZ1
                << ", E = " << Ep[facei] << endl;*/
        }
    }
}


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
        + Erosion
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
