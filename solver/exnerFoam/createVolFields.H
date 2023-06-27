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

// motion solver
/*autoPtr<motionSolver> mSolver;

if (isMeshMoving)
{
    mSolver.reset(motionSolver::New(mesh));
    }*/

// Volume-to surface mapping object
const volSurfaceMapping vsm(aMesh);

const dimensionedVector g(dimensionSet(0, 1, -2, 0, 0, 0, 0), vector(0, 0, -9.81));

volScalarField Zbvf
(
    IOobject
    (
        "Zbvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);

volVectorField Qbvf
(
    IOobject
    (
        "Qbvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0, 2, -1, 0, 0, 0, 0), Zero)
);

vsm.mapToVolume(Zb, Zbvf.boundaryFieldRef());
Zbvf.write();
vsm.mapToVolume(Qb, Qbvf.boundaryFieldRef());
Qbvf.write();
