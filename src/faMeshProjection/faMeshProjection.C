/*---------------------------------------------------------------------------* \
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
    deleteDemandDrivenData(areaCentresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
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
    nInternalEdges_(aMesh.nInternalEdges()),
    nFaces_(aMesh.nFaces()),
    edgeOwner_(aMesh.edgeOwner()),
    edgeNeighbour_(aMesh.edgeNeighbour()),
    projectNormal_(projectNormal/mag(projectNormal)),
    areaCentresPtr_(nullptr),
    edgeCentresPtr_(nullptr)
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

label Foam::faMeshProjection::nInternalEdges() const
{
    return nInternalEdges_;
}

label Foam::faMeshProjection::nFaces() const
{
    return nFaces_;
}

const labelList Foam::faMeshProjection::edgeOwner() const
{
    return edgeOwner_;
}

const labelList Foam::faMeshProjection::edgeNeighbour() const
{
    return edgeNeighbour_;
}

const vectorField& Foam::faMeshProjection::Le() const
{
    if (LePtr_==nullptr)
    {
        calcLe();
    }

    return *LePtr_;
}

const scalarField& Foam::faMeshProjection::magLe() const
{
    if (magLePtr_==nullptr)
    {
        calcMagLe();
    }

    return *magLePtr_;
}

const vectorField& Foam::faMeshProjection::areaCentres() const
{
    if (areaCentresPtr_==nullptr)
    {
        calcAreaCentres();
    }

    return *areaCentresPtr_;
}

const vectorField& Foam::faMeshProjection::edgeCentres() const
{
    if (edgeCentresPtr_==nullptr)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}

const vectorField& Foam::faMeshProjection::pointCoords() const
{
    if (pointCoordsPtr_==nullptr)
    {
        calcPointCoords();
    }

    return *pointCoordsPtr_;
}


const scalarField& Foam::faMeshProjection::S() const
{
    if (SPtr_==nullptr)
    {
        calcS();
    }

    return *SPtr_;
}

vectorField Foam::faMeshProjection::project
(
    vectorField field
) const
{
    return field - (field & projectNormal_) * projectNormal_;
}

vector Foam::faMeshProjection::project
(
    vector vec
) const
{
    return vec - (vec & projectNormal_) * projectNormal_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
