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

#include "readGravitationalAcceleration.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// fluid and sediment properties
dimensionedScalar nuF("nu", transportProperties);
dimensionedScalar rhoS("rhoS", dimDensity, transportProperties);
dimensionedScalar dS("dS", dimLength, transportProperties);
dimensionedScalar Cmax("Cmax", dimless, transportProperties);
dimensionedScalar Sc("Sc", dimless, 1);  // Schmidt number default value
if (transportProperties.found("Sc"))
{
    Sc.value() = readScalar(transportProperties.lookup("Sc"));
}

Info<< "Reading settlingProperties\n" << endl;

IOdictionary settlingProperties
(
    IOobject
    (
        "settlingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// water density, should be changed and read in case
dimensionedScalar rhoF(dimDensity, 1000);

// Initialize Settling Model
Info << "Initialize settling Model\n" << endl;
SettlingModel settlingModel(settlingProperties);

Info<< "Reading field settling velocity Ws\n" << endl;
volVectorField Ws
(
    IOobject
    (
        "Ws",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initialize sediment settling velocity\n" << endl;

Ws = settlingModel.Ufall(C, Cmax, dS, rhoS, rhoF, nuF, mag(g))*g/mag(g);
Ws.correctBoundaryConditions();

Info << "Reading/calculating face flux field phi\n" << endl;

#include "createPhi.H"

Info<< "Reading/calculating face flux field phi\n" << endl;

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
    fvc::interpolate(U+Ws) & mesh.Sf()
);

phip = fvc::flux(U+Ws);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

volVectorField shieldsVf
(
    IOobject
    (
        "shieldsvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, Zero)
);

dimensionedScalar critShields(dimless, 0.047);

volVectorField qbVf
(
    IOobject
    (
        "qbVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

volScalarField thetaVf
(
    IOobject
    (
        "thetaVf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"
