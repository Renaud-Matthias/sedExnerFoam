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

Class
    Foam::faMeshProjection

Description
    Manage the presence of a sediment bed in the domain

SourceFiles
    faMeshProjection.C

\*---------------------------------------------------------------------------*/

#ifndef faMeshProjection_H
#define faMeshProjection_H

#include "faCFD.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMeshProjection
{

private:

    // Primitive data

    const faMesh& mesh_;

    //- Total number of points
    const label nPoints_;

    //- Total number of edges
    const label nEdges_;

    //- Number of internal edges
    const label nInternalEdges_;

    //- Total number of faces
    const label nFaces_;

    //- Edge owner adressing
    const labelList edgeOwner_;

    //- Edge neighbour adressing
    const labelList edgeNeighbour_;

    // Geometric data

    //- Unit vector normal to projected faMesh
    const vector projectNormal_;

    //- Edge length vector oriented from owner to neighbour face
    mutable vectorField* LePtr_;
    
    //- Edge length magnitude
    mutable scalarField* magLePtr_;

    //- Face centres
    mutable vectorField* areaCentresPtr_;

    //- Edge centres
    mutable vectorField* edgeCentresPtr_;

    // Private member functions for demand driven data

    //- Calculate edge length vectors
    void calcLe() const;

    //- Calculate face centres coordinates
    void calcAreaCentres() const;

    //- Calculate edge centres coordinates
    void calcEdgeCentres() const;

    //- Clear geometry
    void clearGeom() const;
    

public:

    // Constructors

    //- Construct from finite area mesh
    faMeshProjection(faMesh& mesh, vector projector);

    // Destructors
    
    ~faMeshProjection();

    // Member functions

    //- Return const reference to finite area mesh
    const faMesh& mesh() const;

    //- Return number of points
    label nPoints() const;

    //- Return number of edges, internal + boundary
    label nEdges() const;

    //- Return number of internal edges
    label nInternalEdges() const;

    //- Return number of faces
    label nFaces() const;

    //- Edge owner adressing
    // ith element is label of face owning edge i
    const labelList edgeOwner() const;

    //- Edge neighbour adressing
    // ith element is label of face neighbour to  edge i
    const labelList edgeNeighbour() const;

    // Demand driven data

    //- Return edge length vectors oriented from owner to neighbour faces
    const vectorField& Le();
    
    //- Return edge length magnitude
    //const scalarField& magLe();

    //- Return projected face centers coordinates
    const vectorField& areaCentres();

    //- Return projected edge centers coordinates
    const vectorField& edgeCentres();

    //- Return edge normals vectors, oriented from owner to neighbour
    //const vectorField& edgeNormals()

    //- Return projected faces area
    //const scalarField& S();
};

}  // End namespace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
