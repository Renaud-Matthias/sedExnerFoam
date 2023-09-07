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
                faBed.aMesh.ref(),
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
                faBed.aMesh.ref(),
                dimensionedVector(dimVelocity * dimLength, Zero)
            )
        );

    faBed.zb.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "zb",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh.ref(),
                dimensionedVector(dimLength, Zero)
            )
        );

    // initiate zb value with fa mesh bottom boundary faces coordinates
    forAll(faBed.aMesh.ref().areaCentres(), facei)
    {
        faBed.zb.ref()[facei] = faBed.aMesh.ref().areaCentres()[facei];
    }

    faBed.deltaH.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "deltaH",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh.ref(),
                dimensionedVector(dimLength, Zero)
            )
        );


    faBed.phib.reset
        (
            new edgeScalarField
            (
                IOobject
                (
                    "phib",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh.ref(),
                dimensionedScalar(dimensionSet(0, 3, -1, 0, 0, 0, 0), Zero)
            )
        );

    faBed.thetaBed.reset
        (
            new areaScalarField
            (
                IOobject
                (
                    "thetaBed",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                faBed.aMesh.ref(),
                dimensionedScalar(dimless, Zero)
            )
        );
}
