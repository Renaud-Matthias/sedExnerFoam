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

    //- Total number of faces
    const label nFaces_;

    // Geometric data

    //- Unit vector normal to projected faMesh
    const vector projectNormal_;

    //- Face centres
    mutable vectorField* faceCentresPtr_;

    // Private member functionsfor demand driven data

    //- Calculate face centres coordinates
    void calcFaceCentres() const;

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

    //- Return number of faces
    label nFaces() const;

    // Demand driven data

    //-
    const vectorField& faceCentres();

    //- Return projected faces area
    //const scalarField& S();
};
    
}  // End namespace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
