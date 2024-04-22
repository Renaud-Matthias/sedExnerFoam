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

// points displacements values for mesh motion
autoPtr<vectorField> dispValsPtr;

// horizontal projection of faMesh
autoPtr<faMeshProjection> aMeshProjPtr;
// fields for bedload and morphology
autoPtr<areaVectorField> shieldsPtr;  // shields number
autoPtr<areaScalarField> critShieldsPtr;  // critical shields number
autoPtr<areaVectorField> qbPtr;  // bedload
autoPtr<areaScalarField> dHptr;  // bed elevation increment on faces
autoPtr<scalarField> dHpointsPtr; // bed elevation increment on points
autoPtr<areaScalarField> DEptr;  // deposition/erosion source term
autoPtr<edgeScalarField> phiqbPtr;  // bedload flux through edges
autoPtr<areaScalarField> dVsedPtr;  // volumic variation, face displacements
autoPtr<surfaceScalarField> bedPhiPtr;  // mesh flux on bed patch
// bedloadModel
autoPtr<bedloadModel> bedloadModelPtr;  // bedloadModel

// repose angle of granular material (in rad)
dimensionedScalar betaRep(dimless, 0.5585);
betaRep.value() = degToRad
    (
        bedloadProperties.lookupOrDefault<scalar>("reposeAngle", 32)
    );

// morphological acceleration factor, default is 1
scalar morphoAccFactor = 1;
morphoAccFactor =
    bedloadProperties.lookupOrDefault<scalar>("morphoAccFactor", 1);
Info << "morphological acceleration factor is set to: "
<< morphoAccFactor << endl;

// coefficient for avalanche additional bedload
dimensionedScalar Qav(dimVelocity * dimLength, 0.01);
Qav.value() = bedloadProperties.lookupOrDefault<scalar>("Qav", 0.01);

// create fa fields if morphology on
if (bed.exist())
{
    Info << "bedload is on" << endl
        << "Instantiate faMeshProjection" << endl;

    aMeshProjPtr.reset
    (
        new faMeshProjection(bed.aMesh(), (g / mag(g)).value())
    );
    
    Info << "Instantiate finite area fields" << endl;

    dispValsPtr.reset(new vectorField(bed.aMesh().nPoints()));
    vectorField& dispVals = dispValsPtr.ref();

    dHpointsPtr.reset(new scalarField(bed.aMesh().nPoints()));

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

    bedPhiPtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "bedPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimensionSet(0, 3, -1, 0, 0, 0, 0), Zero)
        )
    );

    dVsedPtr.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "dVsed",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimensionSet(0, 3, -1, 0, 0, 0, 0), Zero)
        )
    );

    DEptr.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "DE",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimVelocity, Zero)
        )
    );

    if (bedloadProperties.found("bedloadModel"))
    {
        dictionary bedloadModelDict =
            bedloadProperties.subDict("bedloadModel");

        bedloadModelPtr.reset
            (
                new bedloadModel
                (
                    bedloadModelDict
                )
            );
    }
    else
    {
        bedloadModelPtr.reset(new bedloadModel());
    }

    const bedloadModel& bedloadMod = bedloadModelPtr.ref();

    Info << "\nbedload Model: " << endl
        << bedloadMod << endl;

    critShieldsPtr.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "critShields",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            bedloadMod.criticalShields0()
        )
    );
}
