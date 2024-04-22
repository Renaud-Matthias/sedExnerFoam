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

// Access to variables relative to bedload
areaVectorField& shields = shieldsPtr.ref();
areaVectorField& qb = qbPtr.ref();

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

dimensionedScalar smallVal(dimless, SMALL);
// unit vector gravity oriented
vector eg = (g/mag(g)).value();

// Access bedload model and critical Shields number
const bedloadModel& bedloadMod = bedloadModelPtr.ref();
areaScalarField& critShields = critShieldsPtr.ref();
// base critical Shields number
const dimensionedScalar& critShields0 = bedloadMod.criticalShields0();

// shear stress direction
vectorField stressDir = shields / (Foam::mag(shields) + SMALL);
// slope direction
const vectorField& nFaces =
    bed.aMesh().faceAreaNormals().internalField();
vectorField slopeDir = (nFaces & eg) * eg - nFaces;
slopeDir /= (Foam::mag(slopeDir) + SMALL);
forAll(critShields, facei)
{
    scalar beta = Foam::acos(nFaces[facei] & eg);  // bed slope angle
    // angle between shear stress and slope direction
    scalar alpha = Foam::acos(slopeDir[facei] & stressDir[facei]);
    // slope effect on critical Shields number
    scalar betap = Foam::min(beta, betaRep.value());
    critShields[facei] = critShields0.value()
      * (
          Foam::cos(betap)
        * Foam::sqrt
          (
              1 - Foam::pow
              (
                  Foam::sin(alpha)
                * Foam::tan(betap)
                * (1/Foam::tan(betaRep.value())),
                  2
              )
          )
          - Foam::cos(alpha) * Foam::sin(betap) / Foam::tan(betaRep.value())
      );
}

// compute bedload from Shields number
qb = bedloadMod.qb(shields, critShields, rhoS, rhoF, mag(g), dS);

// additional bedload due to avalanche
if (bed.isAvalanche())
{
    #include "sedimAvalanche.H"
}

qb.correctBoundaryConditions();
