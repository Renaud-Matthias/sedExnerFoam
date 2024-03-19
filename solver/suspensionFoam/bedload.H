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

// Access to bed variables
areaVectorField& shields = shieldsPtr.ref();
areaVectorField& qb = qbPtr.ref();
areaScalarField& dH = dHptr.ref();
edgeScalarField& phiqb = phiqbPtr.ref();

vector gUnit = (g / mag(g)).value();  // unit vector gravity oriented

// compute wallShearStress on each patch belonging to sediment bed
for (const label patchID : bed.bedPatchesID())
{
    // compute wall shear stress on patch
    vectorField& ssp = shieldsVf.boundaryFieldRef()[patchID];
    const vectorField& Sfp = mesh.Sf().boundaryField()[patchID];
    const scalarField& magSfp =
        mesh.magSf().boundaryField()[patchID];
    volSymmTensorField Reff =
        rhoF * turbulence->nuEff() * dev(twoSymm(fvc::grad(U)));
    const symmTensorField& Reffp = Reff.boundaryField()[patchID];
    
    // compute shields number
    ssp = ((-Sfp / magSfp) & Reffp)
        / (mag(g).value() * dS.value()
        * (rhoS.value() - rhoF.value()));
}

// map volumic shields number to area shields number
shields.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface<vector>
    (
        shieldsVf.boundaryFieldRef()
    );

// compute bedload with Meyer Peter Muller law

//dimensionedScalar mpmE(dimless, 8);
dimensionedScalar smallVal(dimless, SMALL);

// compute bedload from Shields number
const bedloadModel& bedloadMod = bedloadModelPtr.ref();

qb = bedloadMod.qb(shields, rhoS, rhoF, mag(g), dS);

/*forAll(bed.aMesh().areaCentres(), facei)
{
    dimensionedScalar shieldsFacei(dimless, shields[facei]);
    qb[facei] = (shields[facei] / (smallVal + Foam::mag(shields[facei])))
        * bedloadMod.qb
        (
            mag(shieldsFacei),
            rhoS,
            rhoF,
            Foam::mag(g),
            dS
        ).value();
}*/

/*qb =
    mpmE * (shields / (smallVal + mag(shields)))
    * Foam::sqrt((((rhoS / rhoF) - 1) * mag(g) * Foam::pow(dS, 3)))
    * Foam::pow(
        pos(mag(shields) - critShields)
    * (mag(shields) - critShields), 1.5);*/

qb.correctBoundaryConditions();

if (bed.isAvalanche())
{
    #include "sedimAvalanche.H"
}

// compute qb flux through fa mesh edges
phiqb =
    linearEdgeInterpolate(qb) & bed.aMesh().Le();

// solve exner explicitly (euler explicit)
// bedload contribution
dH = -runTime.time().deltaT()
    * morphoAccFactor
    * (1/CsMax)
    * (1/(bed.aMesh().faceAreaNormals() & (gUnit)))  // slope correction
    * fac::edgeIntegrate(phiqb);

if (switchSuspension=="on")
{
    #include "erosionDeposition.H"
}

if (bed.meshMotion())
{
    #include "meshMove.H"
}
