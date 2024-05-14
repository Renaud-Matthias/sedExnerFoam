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

// Instantiate sedimentBed object
Info << "Instantiate sedimentBed" << endl;
sedimentBed bed(bedloadProperties, mesh, g);

// Check bedload state (on/off)
// bedload is off if no sedimentBed exists
word switchBedload("off");
if (bed.exist())
{
    switchBedload = word(bedloadProperties.lookup("bedload"));
}
if (switchBedload!="on" and switchBedload!="off")
{
    FatalError << switchBedload
        << " is not a valid entry for bedload"
        << ", valid options are: on off" << endl;
    Info << abort(FatalError) << endl;
}

Info << "bedload is " << switchBedload << endl;

// points displacements values for mesh motion
autoPtr<vectorField> dispValsPtr;

// horizontal projection of faMesh
autoPtr<faMeshProjection> aMeshProjPtr;
// fields, variables for bedload and morphology
autoPtr<areaVectorField> shieldsPtr;  // shields number
autoPtr<areaScalarField> critShieldsPtr;  // critical shields number
autoPtr<areaVectorField> qbPtr;  // bedload
autoPtr<areaScalarField> dHptr;  // bed elevation increment on faces
autoPtr<scalarField> dHpointsPtr; // bed elevation increment on points
autoPtr<areaScalarField> DEptr;  // deposition/erosion source term
autoPtr<edgeScalarField> phiqbPtr;  // bedload flux through edges
autoPtr<areaScalarField> dVsedPtr;  // volumic variation, face displacements
autoPtr<surfaceScalarField> bedPhiPtr;  // mesh flux on bed patch
// models relative to bedload and morphology
autoPtr<bedloadModel> bedloadModelPtr;  // bedloadModel
// critical Shields
autoPtr<criticalShieldsModels::criticalShieldsModel>
critShieldsModelPtr;

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
    Info << "Instantiate faMeshProjection" << endl;

    aMeshProjPtr.reset
    (
        new faMeshProjection(bed.aMesh(), (g / mag(g)).value())
    );

    // vertices displacements for mesh motion
    dispValsPtr.reset(new vectorField(bed.aMesh().nPoints()));
    vectorField& dispVals = dispValsPtr.ref();
    // bed vertices vertical displacements, Exner equation
    dHpointsPtr.reset(new scalarField(bed.aMesh().nPoints()));
    for (label pointi = 0; pointi < bed.aMesh().nPoints(); pointi++)
    {
        dispVals[pointi] = vector(0, 0, 0);
    }


    // model for critical Shields number
    if (bedloadProperties.found("criticalShieldsModel"))
    {
        const dictionary critShieldsDict =
            bedloadProperties.subDict("criticalShieldsModel");

        critShieldsModelPtr.reset
        (
            criticalShieldsModels::criticalShieldsModel::New
            (
                critShieldsDict
            )
        );
    }


    // model for bedload
    if (bedloadProperties.found("bedloadModel"))
    {
        const dictionary bedloadModelDict =
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

    Info << "Instantiate finite area fields" << endl;
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

    // add a turbulent diffusivity for suspended particles near wall
    // initialize a wallDist function for ditance to bed patches
    if (switchSuspension == "on")
    {
        labelHashSet bedPatchSet;

        for (const label patchID : bed.bedPatchesID())
        {
            bedPatchSet.insert(patchID);
        }

        bedDistPtr.reset(new wallDist(mesh, bedPatchSet, "wall"));

        wallDist& bedDist = bedDistPtr.ref();

        Info << "bed patches ID: " << bedDist.patchIDs() << endl;
        Info << "bed dist y: " << bedDist.y() << endl;
    }
}
