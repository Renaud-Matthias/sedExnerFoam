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
sedbedManager bed(bedloadProperties, mesh, g);

autoPtr<vectorField> dispValsPtr;

// fields for bedload and morphology
autoPtr<areaVectorField> shieldsPtr;
autoPtr<areaVectorField> qbPtr;
autoPtr<areaScalarField> dHptr;
autoPtr<edgeScalarField> phiqbPtr;

// create fa fields if morphology on
if (bed.exist())
{
    Info << "Instantiate finite area fields" << endl;

    dispValsPtr.reset(new vectorField(bed.aMesh().nPoints()));
    
    vectorField& dispVals = dispValsPtr.ref();

    for (label pointi = 0; pointi < bed.aMesh().nPoints(); pointi++)
    {
        dispVals[pointi] = vector(0, 0, 0);
    }
    
    shieldsPtr.reset
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
                bed.aMesh(),
                dimensionedVector(dimless, Zero)
            )
        );

    qbPtr.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "qb",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                bed.aMesh()
            )
        );

    dHptr.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "deltaH",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimLength, Zero)
        )
    );

    phiqbPtr.reset
    (
        new edgeScalarField
        (
            IOobject
            (
                "phiqb",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimensionSet(0, 3, -1, 0, 0, 0, 0), Zero)
        )
    );
}
