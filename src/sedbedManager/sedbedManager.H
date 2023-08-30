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
    Foam::sedbedManager

Description
    Manage the presence of a sediment bed in the domain

SourceFiles
    sedbedManager.C

\*---------------------------------------------------------------------------*/

#ifndef sedbedManager_H
#define sedbedManager_H

#include "fvCFD.H"
#include "dictionary.H"
#include "dynamicFvMesh.H"
#include "meshTools.H"
#include "faCFD.H"
#include "Time.H"
#include "gravityMeshObject.H"
#include "pointMesh.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

class sedbedManager
{
public:

    autoPtr<faMesh> aMesh_;
    autoPtr<areaVectorField> shields;  // shields number
    autoPtr<areaVectorField> qb;  // bedload transport
    autoPtr<areaVectorField> zb;  // bed elevation
    autoPtr<areaVectorField> deltaH;  // bed level increment
    autoPtr<edgeScalarField> phib;  // qb flux through fa mesh edges
    autoPtr<pointVectorField> pointDisp; // point displacement for mesh motion
    autoPtr<volSurfaceMapping> vsm;  // to read faFields with paraview

     //- Constructor
    sedbedManager
    (
        dictionary& dict,
        fvMesh& mesh,
        const meshObjects::gravity& g
    );


    //- Destructor
    virtual ~sedbedManager();

    //- Public member functions

    // Return copy of patches ID labelList
    labelList bedPatchesID();

    // Check if bed exists
    void checkBedExistence();

    // Existence of sediment bed
    bool exist() const;

    // Mesh motion activated or not
    bool meshMotion() const;
    
protected:

    bool bedExist_;
    bool meshMove_;
    const dictionary& dict_;
    fvMesh& mesh_;
    const meshObjects::gravity& g_;  // gravitational acceleration
    // List of patches names corresponding to a sediment bed
    List<word> bedPatchesNames_;
    labelList bedPatchesID_;

    //- Protected member functions

    // Check if patch are well oriented
    void checkFaMeshOrientation() const;

    // get sediment bed patches ID
    void getPatchesID();
    
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
