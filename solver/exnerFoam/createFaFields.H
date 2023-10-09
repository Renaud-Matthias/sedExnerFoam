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

// Create Finite Area mesh
faMesh aMesh(mesh);
label Psize = aMesh.points().size();

//mesh motion solver
//autoPtr<motionSolver> motionPtr;

Info << "aMesh n points : " << Psize << endl;

// Create Point mesh
pointMesh pMesh(mesh);

Info<< "Reading field Zbf" << endl;
areaScalarField Zb
(
    IOobject
    (
        "Zb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

areaVectorField deltaH
(
    IOobject
    (
        "deltaH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedVector(dimLength, Zero)
);

edgeScalarField Zbedges
(
    IOobject
    (
        "Zbedges",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Zb)
);

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

#include "getSetUp.H"

// moving mesh parameters
bool isMeshMoving(false);
word movingMesh(setupProperties.lookup("movingMesh"));

if (movingMesh=="on")
{
    isMeshMoving = true;
    Info << "mesh motion enable" << endl;
}
else if (movingMesh!="off")
{
    FatalError << "wrong keyword movingMesh in setupProperties" << endl
        << "options are : yes  no" << endl;
    Info << abort(FatalError) << endl;
}

areaVectorField Qb
(
    IOobject
    (
        "Qb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

areaVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedVector(dimVelocity, vector(0, 0, 0))
);

dimensionedScalar dLength(dimLength, SMALL);
edgeScalarField phiqb
(
    IOobject
    (
        "phiqb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Qb) & aMesh.Le()
);

edgeScalarField phiub
(
    IOobject
    (
        "phiub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Ub) & aMesh.Le()
);

areaScalarField Da
(
    IOobject
    (
        "Da",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

areaScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimless, Zero)
);

areaScalarField rhsExnerEq
(
    IOobject
    (
        "rhsExnerEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimensionSet(0, 1, -1, 0, 0, 0, 0), Zero)
);

/* ---- fields for filtering exner equation ---- */

word filter(setupProperties.lookup("filter"));
word filterType;

if (filter=="on")
{
    filterType = word(setupProperties.lookup("filterType"));
}

// number of loop for bed shape corrector
scalar alphaFilter = readScalar(setupProperties.lookup("alphaFilter"));
dimensionedScalar diffFilter(dimensionSet(0, 2, 0, 0, 0, 0, 0), 0);
const label Ncorr(readLabel(setupProperties.lookup("Ncorrection")));

areaScalarField Zbpred
(
    IOobject
    (
        "Zbpred",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimLength, Zero),
    "zeroGradient"
);

areaScalarField dzCorr
(
    IOobject
    (
        "dzCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimLength, Zero),
    "zeroGradient"
);
