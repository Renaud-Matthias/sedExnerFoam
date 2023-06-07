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

Info<< "Reading bedloadProperties\n" << endl;

IOdictionary bedloadProperties
(
    IOobject
    (
        "bedloadProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Instantiate sedbedManager object
Info << "Instantiate sedbedManager" << endl;
sedbedManager faBed(bedloadProperties, mesh, g);

//Info << "check sediment bed existence" << endl;

if (faBed.exist())
{
    Info << "Instantiate finite area fields" << endl;
    
    faBed.shields.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "shields",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh_.ref(),
                dimensionedVector(dimless, Zero)
            )
        );

    faBed.qb.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "qb",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh_.ref(),
                dimensionedVector(dimVelocity * dimLength, Zero)
            )
        );
}
