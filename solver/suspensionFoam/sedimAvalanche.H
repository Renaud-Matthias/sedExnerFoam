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

// compute bed slope angle
const vectorField& nFaces =
    bed.aMesh().faceAreaNormals().internalField();

//bed.angleBedSlope().primitiveFieldRef() = Foam::acos(nFaces & gUnit);

// slope direction
vectorField slopeDir = (nFaces & gUnit) * gUnit - nFaces;
slopeDir /= (mag(slopeDir) + SMALL);

// reposeAngle
//dimensionedScalar betaRep(dimless, degToRad(32));
//dimensionedScalar Qav(dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.01);

// Vinent law (2019)
forAll(bed.aMesh().areaCentres(), facei)
{
    scalar beta = Foam::acos(nFaces[facei] & gUnit);
    vector qav = Qav.value() * slopeDir[facei]
        * pos(beta - betaRep.value())
        * (
            Foam::tanh(Foam::tan(beta))
            - Foam::tanh(Foam::tan(betaRep.value())));
    qb[facei] += qav;
}
