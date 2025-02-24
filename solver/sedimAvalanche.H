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

// unit vector gravity oriented
vector eg = (g/mag(g)).value();

// compute bed slope angle
const vectorField& nFaces =
    bed.aMesh().faceAreaNormals().internalField();

//bed.angleBedSlope().primitiveFieldRef() = Foam::acos(nFaces & gUnit);

// bed faces angle
scalarField slopeAngle = Foam::acos(nFaces & eg);

// slope direction
vectorField slopeDir = eg - nFaces * (eg & nFaces);
slopeDir /= (Foam::mag(slopeDir) + SMALL);

// Vinent law (2019)
vectorField& qbIntern = qb.primitiveFieldRef();
qbIntern += bedloadModelPtr->qbAvalanche(slopeAngle, slopeDir, betaRep);
/*forAll(bed.aMesh().areaCentres(), facei)
{
    scalar beta = Foam::acos(nFaces[facei] & eg);  // bed slope angle
    vector qav = Qav.value() * slopeDir[facei]
        * pos(beta - betaRep.value())
        * (
            Foam::tanh(Foam::tan(beta))
            - Foam::tanh(Foam::tan(betaRep.value())));
    qb[facei] += qav;
    }*/
