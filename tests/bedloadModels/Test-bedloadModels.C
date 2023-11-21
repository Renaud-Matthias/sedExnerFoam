/*---------------------------------------------------------------------------* \
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

#include "fvCFD.H"
#include "bedloadModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "Run tests for bedload Model" << endl;

    bedloadModel mpm("Meyer-Peter");

    bedloadModel nielsen("Nielsen");

    scalar rhoS = 2500;
    scalar rhoF = 1000;
    scalar g = 9.81;
    scalar dS = 0.001;
    scalarField shieldsField(5);
    forAll(shieldsField, vali)
    {
        scalar val = 0.02 * vali;
        shieldsField[vali] = val;
    }
    

    Info << mpm;
    Info << "for shields = " << shieldsField
        << ", qb = " << mpm.qb(shieldsField, rhoS, rhoF, g, dS) << endl;

    Info << nielsen;
    Info << "for shields = " << shieldsField
        << ", qb = " << nielsen.qb(shieldsField, rhoS, rhoF, g, dS) << endl;

    return 0;
}

// ************************************************************************* //
