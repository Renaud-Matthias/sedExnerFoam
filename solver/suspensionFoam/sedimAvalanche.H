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

// slope direction
vectorField slopeDir = (nFaces & gUnit) * gUnit - nFaces;
slopeDir /= (mag(slopeDir) + SMALL);

// reposeAngle
dimensionedScalar thetaRep(dimless, degToRad(32));
dimensionedScalar Qav(dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.01);
// Vinent law (2019)
forAll(faBed.aMesh.ref().areaCentres(), i)
{
    scalar theta = faBed.thetaBed.ref()[i];
    vector qav = Qav.value() * slopeDir[i]
        * pos(theta - thetaRep.value())
        * (
            Foam::tanh(Foam::tan(theta))
            - Foam::tanh(Foam::tan(thetaRep.value())));
    faBed.qb.ref()[i] += qav;
}
