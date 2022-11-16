Info<< "Reading field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading sedimentProperties\n" << endl;

IOdictionary sedimentProperties
(
    IOobject
    (
        "sedimentProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoS("rhoS", dimDensity, sedimentProperties);
dimensionedScalar dS("dS", dimLength, sedimentProperties);

settlingModel UfallModel(sedimentProperties);
/*settlingModel Ufallmodel(sedimentProperties);
autoPtr<settlingModel> pp = settlingModel::New
(
    sedimentProperties
);*/

Info<< "Compute sediment settling velocity Ufall\n" << endl;

dimensionedVector Ufall("Ufall", dimVelocity, vector (0, -0.05, 0));

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading/calculating face flux field phi and phip\n" << endl;

#include "createPhi.H"

surfaceScalarField phip
(
    IOobject
    (
        "phip",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U+Ufall)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
