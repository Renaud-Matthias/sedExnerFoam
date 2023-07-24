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

areaScalarField deltaH
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
    dimensionedScalar(dimLength, Zero)
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

if (movingMesh=="yes")
{
    isMeshMoving = true;
}
else if (movingMesh!="no")
{
    FatalError << "wrong keyword movingMesh in transportProperties" << endl
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

edgeScalarField phib
(
    IOobject
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Qb) & aMesh.Le()
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
    dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), Zero)
);

/* ---- fields for filtering exner equation ---- */

// number of loop for bed shape corrector
const label Nouter(readLabel(setupProperties.lookup("Nouter")));
const label Ninner(readLabel(setupProperties.lookup("Ninner")));

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

scalar alphaFilter = readScalar(setupProperties.lookup("alphaFilter"));

areaScalarField diffFilter
(
    IOobject
    (
        "diffFilter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedScalar(dimArea, Zero),
    "zeroGradient"
);
