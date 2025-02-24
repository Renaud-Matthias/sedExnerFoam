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

// Access to saturated bedload and bedload fields
areaVectorField& qsat = qsatPtr.ref();
areaVectorField& qb = qbPtr.ref();

dimensionedScalar smallVal(dimless, SMALL);
// unit vector gravity oriented
vector eg = (g/mag(g)).value();

// Access bedload model and critical Shields number
const bedloadModels::bedloadModel& bedload = bedloadModelPtr.ref();
areaScalarField& critShields = critShieldsPtr.ref();

if (not readCritShields)
{
    // shear stress direction
    vectorField stressDir = shields / (Foam::mag(shields) + SMALL);
    // slope direction
    const vectorField& nFaces =
        bed.aMesh().faceAreaNormals().internalField();
    // bed faces angle
    scalarField slopeAngle = Foam::acos(nFaces & eg);
    // slope direction
    vectorField slopeDir = eg - nFaces * (eg & nFaces);
    slopeDir /= (Foam::mag(slopeDir) + SMALL);
    //vectorField slopeDir = (nFaces & eg) * eg - nFaces;
    //slopeDir /= (Foam::mag(slopeDir) + SMALL);
    // Apply slope correction on critical Shields number
    critShieldsModelPtr->slopeCorrection(
        critShields, stressDir, slopeAngle, slopeDir, betaRep);
}

// compute bedload from Shields number
qsat = bedload.qb(shields, critShields, rhoS, rhoF, mag(g), dS);

// apply saturation time to compute bedload
const scalar& deltaT = runTime.time().deltaT().value();
if (Tsat)
{
    Info << "Tsat: " << Tsat << endl;
    if (Tsat->value() > deltaT)
    {
        qb += (deltaT / Tsat->value()) * (qsat - qb);
    }
}
else
{
    qb = qsat;
}

// additional bedload due to avalanche
if (bedloadModelPtr->avalanche())
{
    #include "sedimAvalanche.H"
}

qb.correctBoundaryConditions();
