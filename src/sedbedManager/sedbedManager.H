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

#include "dictionary.H"
#include "dynamicFvMesh.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

class sedbedManager
{
protected:

    bool bedExist_;
    const dictionary& dict_;
    fvMesh& mesh_;
    const meshObjects::gravity& g_;  // gravitational acceleration
    // List of patches names corresponding to a sediment bed
    List<word> bedPatchesNames_;
    List<label> bedPatchesID_;
    
public:

    //- Constructor
    sedbedManager
    (
        dictionary& dict,
        fvMesh& mesh,
        const meshObjects::gravity& g
    );


    //- Destructor
    virtual ~sedbedManager();

    //- Member functions
    // Existence of sedbed patches
    bool exist() const;

    // Return
    void findBedPatches();

    // Return List of sedbed patches names
    List<word> getPatchesNames() const;

    // Check patch existence and return patch ID
    label checkPatchExistence(word patchName) const;
    
    // Check if patch are well oriented
    void checkPatchOrientation(const polyPatch& aPatch) const;
    
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
