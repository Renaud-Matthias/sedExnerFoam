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

// points displacements values for mesh motion
autoPtr<vectorField> dispValsPtr;

// horizontal projection of faMesh
autoPtr<faMeshProjection> aMeshProjPtr;
// fields, variables for bedload and morphology
autoPtr<areaVectorField> shieldsPtr;  // shields number
autoPtr<areaScalarField> critShieldsPtr;  // critical shields number
autoPtr<areaVectorField> qsatPtr;  // saturated bedload
autoPtr<areaVectorField> qbPtr;  // bedload
autoPtr<areaScalarField> dHptr;  // bed elevation increment on faces
autoPtr<scalarField> dHpointsPtr; // bed elevation increment on points
autoPtr<areaScalarField> DEptr;  // deposition/erosion source term
autoPtr<scalarField> phiqbProjPtr;  // projected bedload flux through edges
autoPtr<areaScalarField> dVsedPtr;  // volumic variation, face displacements
autoPtr<surfaceScalarField> bedPhiPtr;  // mesh flux on bed patch
autoPtr<areaVectorField> rigidBedPtr;  // rigid bed coordinates
// width of sediment bed layer, between the rigid bed and the fluid
autoPtr<areaScalarField> HsedBedPtr;

// models to compute bedload
autoPtr<bedloadModels::bedloadModel> bedloadModelPtr;
// Pointer to critical Shields model
autoPtr<criticalShieldsModels::criticalShieldsModel>
critShieldsModelPtr;
// True if critShields is instantiated from time/finite-area/critShields
bool readCritShields = false;

// saturated length and time
dimensionedScalar* Tsat = nullptr;  // pointer to saturation time
if (bedloadProperties.found("saturation"))
{
    const dictionary& satDict = bedloadProperties.subDict("saturation");
    Tsat = new dimensionedScalar(
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        readScalar(satDict.lookup("Tsat")));
    Info << "saturation time, Tsat = " << Tsat->value() << " s" << endl;
}

scalar betaRepDeg =
    (
        transportProperties.lookupOrDefault<scalar>("reposeAngle", 32)
    );

Info << "granular material repose angle: " << betaRepDeg
    << " degrees" << endl;
if (betaRepDeg <= 0 or betaRepDeg >= 90)
{
    FatalError
        << "repose angle must be in the interval 0 - 90" << endl;
    Info << abort(FatalError) << endl;
}
scalar betaRep = degToRad(betaRepDeg);

// morphological acceleration factor, default is 1
scalar morphoAccFactor = 1;
morphoAccFactor =
    bedloadProperties.lookupOrDefault<scalar>("morphoAccFactor", 1);
Info << "morphological acceleration factor is set to: "
<< morphoAccFactor << endl;

// coefficient for avalanche additional bedload
dimensionedScalar Qav(dimVelocity * dimLength, 0.01);
Qav.value() = bedloadProperties.lookupOrDefault<scalar>("Qav", 0.01);

// list of patch to apply zeroDisplacements condition
List<word> fixedFaPatchNames;

if (bedloadProperties.found("fixedFaPatches"))
{
    fixedFaPatchNames = bedloadProperties.lookup("fixedFaPatches");
}
Info << "fixedFaPatchesNames: " << fixedFaPatchNames << endl;

List<label> fixedFaPatchLabels;
autoPtr<labelList> fixedPatchPointLabelsPtr;

// order for explicit Adams-Bashforth time scheme
// default is 0, euler explicit
label ABorder = 0;
word qbInterpScheme = "linear";
bool qbUpwind = false;

// create fa fields if morphology on
if (bed.exist())
{
    const faBoundaryMesh& faBounds = bed.aMesh().boundary();
    const uindirectPrimitivePatch& bp = bed.aMesh().patch();
    label nFixedPatches = fixedFaPatchNames.size();
    Info << "\nnumber of fixed faPatches: " << nFixedPatches << endl;
    for (label bci=0; bci<nFixedPatches; bci++)
    {
        const word& faPatchName = fixedFaPatchNames[bci];
        label faPatchID = faBounds.findPatchID(fixedFaPatchNames[bci]);
        if (faPatchID==-1)
        {
            FatalError << "in bedloadProperties,"
                << "wrong entry for fixedFaPatches " << faPatchName
                << "\nvalid options are: " << faBounds.names() << endl;
            Info << abort(FatalError) << endl;
        }
        fixedFaPatchLabels.append(faPatchID);
    }
    Info << "fixedFaPatchID: " << fixedFaPatchLabels << endl;

    // create a labelHashSet to store faPatch vertices labels
    labelHashSet verticesSet;
    
    for (label bci=0; bci<fixedFaPatchLabels.size(); bci++)
    {
        label patchID = fixedFaPatchLabels[bci];
        const faPatch& patch = faBounds[patchID];
        const labelList& edgeList = patch.edgeLabels();
        for (const label edgei : edgeList)
        {
            // edge with local addressing
            const edge& e = bp.edges()[edgei];
            verticesSet.insert(e.start());
            verticesSet.insert(e.end());
        }
    }
    
    fixedPatchPointLabelsPtr.reset(new labelList(verticesSet.size()));
    Info << "vertices labelHashSet: " << verticesSet << endl;
    labelList& fixedPatchPointLabels = fixedPatchPointLabelsPtr.ref();

    label icount = 0;
    for (const label pointLabel : verticesSet)
    {
        Info << pointLabel << "; " << endl;
        fixedPatchPointLabels[icount] = pointLabel;
        icount ++;
    }
    
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
    dictionary defaultCritShieldsDict;
    defaultCritShieldsDict.add(keyType("type"), "fixedValue");
    defaultCritShieldsDict.add(keyType("value"), 0.047);
    // if dictionary not present, defaultCritShieldsDict is applied
    // default is fixedValue, value = 0.047
    if (bedloadProperties.found("criticalShieldsModel"))
    {
        const dictionary critShieldsDict =
            bedloadProperties.subDict("criticalShieldsModel");
        word modelType(critShieldsDict.lookup("type"));
        if (modelType != "readFromFile")
        {
            critShieldsModelPtr.reset
                (
                    criticalShieldsModels::criticalShieldsModel::New
                    (
                        critShieldsDict
                    )
                );
            Info << "critical Shields number base value: "
                << critShieldsModelPtr->criticalShields0(Dstar).value()
                << endl;
        }
        else
        {
            readCritShields = true;
            Info << "critical Shields number set by user as "
                << " a non uniform scalar field" << endl;
        }
    }
    else
    {
        critShieldsModelPtr.reset
        (
            criticalShieldsModels::criticalShieldsModel::New
            (
                defaultCritShieldsDict
            )
        );
        Info << "critical Shields number base value: "
            << critShieldsModelPtr->criticalShields0(Dstar).value()
            << endl;
    }

    // model for bedload
    dictionary defaultBedloadDict;
    defaultBedloadDict.add(keyType("type"), "MeyerPeter");
    // if dictionary not present, defaultBedloadDict is applied
    // default is formula from Meyer Peter
    if (bedloadProperties.found("bedloadModel"))
    {
        const dictionary bedloadModelDict =
            bedloadProperties.subDict("bedloadModel");

        bedloadModelPtr.reset
        (
            bedloadModels::bedloadModel::New
            (
                bedloadModelDict
            )
        );
    }
    else
    {
        bedloadModelPtr.reset(
            bedloadModels::bedloadModel::New
            (
                defaultBedloadDict
            )
        );
    }

    // choose time scheme for Exner equation, default is euler explicit
    if (bedloadProperties.found("ABorder"))
    {
        ABorder = readLabel(
            bedloadProperties.lookup("ABorder"));
        if (ABorder!=0 and ABorder!=1)
        {
            FatalError << "in bedloadProperties,"
                << "wrong time scheme entry ABorder " << ABorder
                << "\nvalid options are: 0 1" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    // choose face to edges interpolation scheme for bedload
    qbInterpScheme =
        bedloadProperties.lookupOrDefault<word>(
            "bedloadInterpolation", "linear");
    if (qbInterpScheme != "linear" and qbInterpScheme != "upwind")
    {
        FatalError << "in bedloadProperties,"
            << "wrong entry for bedloadInterpolation " << qbInterpScheme
            << "\nvalid options are: linear upwind" << endl;
        Info << abort(FatalError) << endl;
    }
    if (qbInterpScheme == "upwind")
    {
        qbUpwind = true;
    }
    

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

    qsatPtr.reset
    (
        new areaVectorField
        (
            IOobject
            (
                "qsat",
                runTime.timeName(),
                bed.aMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
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
                bed.aMesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            bed.aMesh()
        )
    );

    phiqbProjPtr.reset
        (
            new scalarField(bed.aMesh().nEdges(), 0)
        );

    dHptr.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "deltaH",
                runTime.timeName(),
                bed.aMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimLength, Zero)
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
                bed.aMesh(),
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
                bed.aMesh(),
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
                bed.aMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bed.aMesh(),
            dimensionedScalar(dimVelocity, Zero)
        )
    );

    if (readCritShields)
    {
        Info << "critical Shields number read from /finite-area/critShields"
            << "\nno criticalShieldsModel used" << endl;

        critShieldsPtr.reset
            (
                new areaScalarField
                (
                    IOobject
                    (
                        "critShields",
                        runTime.timeName(),
                        bed.aMesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    bed.aMesh()
                )
            );
    }
    else
    {
        critShieldsPtr.reset
            (
                new areaScalarField
                (
                    IOobject
                    (
                        "critShields",
                        runTime.timeName(),
                        bed.aMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    bed.aMesh(),
                    critShieldsModelPtr->criticalShields0(Dstar)
                )
            );
    }

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

        //wallDist& bedDist = bedDistPtr.ref();
        //Info << "bed patches ID: " << bedDist.patchIDs() << endl;
        //Info << "bed dist y: " << bedDist.y() << endl;
    }

    if (bed.rigidBed())
    {
        rigidBedPtr.reset
            (
                new areaVectorField
                (
                    IOobject
                    (
                        "rigidBed",
                        runTime.timeName(),
                        bed.aMesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    bed.aMesh()
                )
            );

        HsedBedPtr.reset
            (
                new areaScalarField
                (
                    IOobject
                    (
                        "HsedBed",
                        runTime.timeName(),
                        bed.aMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    bed.aMesh(),
                    dimensionedScalar(dimLength, Zero)
                )
            );

        areaVectorField& rigidBed = rigidBedPtr.ref();
        areaScalarField& HsedBed = HsedBedPtr.ref();

        HsedBed =
            (
                rigidBed
                - bed.aMesh().areaCentres()
            ) & (g/mag(g)).value();

        Info << "width of sediment layer, max: "
            << max(HsedBed).value() << " m, min: "
            << min(HsedBed).value() << " m, mean: "
            << average(HsedBed).value() << endl;
    }
}
