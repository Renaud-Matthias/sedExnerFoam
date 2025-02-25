/*---------------------------------------------------------------------------* \
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    sedExnerFOAM is free software: you can redistribute it and/or modify it
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
    Foam::sedimentBed

Description
    Manage the presence of a sediment bed in the domain

SourceFiles
    sedbedManager.C

\*---------------------------------------------------------------------------*/

#ifndef sedimentbed_H
#define sedimentbed_H

#include "fvCFD.H"
#include "dictionary.H"
#include "dynamicFvMesh.H"
#include "meshTools.H"
#include "faCFD.H"
#include "faMeshProjection.H"
#include "Time.H"
#include "gravityMeshObject.H"
#include "pointMesh.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

class sedimentBed
{
public:
    autoPtr<volSurfaceMapping> vsm;  // to read faFields with paraview
    // variable for exner filterering

     //- Constructor
    sedimentBed
    (
        dictionary& dict,
        fvMesh& mesh,
        const meshObjects::gravity& g
    );


    //- Destructor
    virtual ~sedimentBed();

    // Public member functions

    //- Return const reference to faMesh
    const faMesh& aMesh();

    //- Return const reference to faMeshProjection
    const faMeshProjection& aProjMesh();
    
    //- Avalanche module activated or not
    //bool isAvalanche() const;

    //- Return copy of patches ID labelList
    labelList bedPatchesID();

    //- Existence of sediment bed
    bool exist() const;

    //- Mesh motion activated or not
    bool bedMotion() const;

    //- Presence of a rigid bed below the sediment layer
    bool rigidBed() const;

    //- Return reference to beta, the bed slope angle of faces in rad
    const scalarField& beta() const;

    //- Return reference to slope direction scalarField
    const vectorField& slopeDir() const;

    //- Interpolate vertical bed displacement from faces center to vertices
    // interpolation scheme from Jacobsen 2015
   void interpFaceToVertices
    (
        scalarField& dHfaces,
        scalarField& dHpoints
    );


protected:

    // protected attributes

    const dictionary& dict_;

    bool bedExist_;
    
    bool bedMotion_;

    //bool avalanche_;

    bool rigidBed_;

    fvMesh& mesh_;

    autoPtr<faMesh> aMesh_;

    autoPtr<faMeshProjection> aProjMesh_;

    const meshObjects::gravity& g_;  // gravitational acceleration

    const vector eg_;  // unit vector oriented with g

    labelList bedPatchesID_;
    
    List<word> bedPatchesNames_;

    mutable scalarField* betaPtr_;

    mutable vectorField* slopeDirPtr_;

    // Protected member functions

    //- check if avalanche model is used
    //void checkAvalancheModel_();
    
    //- Check if bed exists
    void checkBedExistence_();

    //- Check if patch are well oriented
    void checkFaMeshOrientation_() const;

    //- Check if bed/mesh motion (on/off)
    void checkBedMotion_();

    //- get sediment bed patches ID
    void getPatchesID();

    //- Compute beta, the bed slope angle, pointed by betaPtr_
    void calcBeta() const;

    //- Compute slope direction unit vector for each face
    void calcSlopeDir() const;
    
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
