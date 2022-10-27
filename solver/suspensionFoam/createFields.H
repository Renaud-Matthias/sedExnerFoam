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

#include "settlingVelocity.H"

Info<< "Reading diffusivity DT\n" << endl;
dimensionedScalar DT("DT", dimViscosity, transportProperties);

Info<< "Reading/calculating face flux field phi and phip\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

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

#include "createFvOptions.H"
