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

#include "faMeshProjection.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMeshProjection::clearGeom() const
{
    deleteDemandDrivenData(faceCentresPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshProjection::faMeshProjection
(
    faMesh& aMesh,
    vector projectNormal
)
:
    mesh_(aMesh),
    nPoints_(aMesh.nPoints()),
    nEdges_(aMesh.nEdges()),
    nFaces_(aMesh.nFaces()),
    projectNormal_(projectNormal/mag(projectNormal)),
    faceCentresPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMeshProjection::~faMeshProjection()
{
    clearGeom();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const faMesh& Foam::faMeshProjection::mesh() const
{
    return mesh_;
}

label Foam::faMeshProjection::nPoints() const
{
    return nPoints_;
}

label Foam::faMeshProjection::nEdges() const
{
    return nEdges_;
}

label Foam::faMeshProjection::nFaces() const
{
    return nFaces_;
}

const vectorField& Foam::faMeshProjection::faceCentres()
{
    if (faceCentresPtr_==nullptr)
    {
        calcFaceCentres();
    }

    return *faceCentresPtr_;
}

const vectorField& Foam::faMeshProjection::edgeCentres()
{
    if (edgeCentresPtr_==nullptr)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
