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

#include "sedbedManager.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedbedManager::sedbedManager
(
    dictionary& dict,
    Foam::fvMesh& mesh,
    const meshObjects::gravity& g
)
:
    bedExist_(false), meshMotion_(false), dict_(dict), mesh_(mesh), g_(g)
{
    checkBedExistence_();
    if (bedExist_)
    {
        aMesh_.reset(new faMesh(mesh_));
        vsm.reset(new volSurfaceMapping(aMesh_.ref()));
        getPatchesID();
        checkFaMeshOrientation_();
        
        Info << "faMesh patches ID : " << bedPatchesID_ << endl;
        Info << "faMesh patches names : " << bedPatchesNames_ << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedbedManager::~sedbedManager()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// public member functions

labelList Foam::sedbedManager::bedPatchesID()
{
    return bedPatchesID_.clone();
}

bool Foam::sedbedManager::exist() const
{
    if (bedExist_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::sedbedManager::meshMotion() const
{
    if (meshMotion_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//- Protected member functions

void Foam::sedbedManager::checkBedExistence_()
{
    if (dict_.found("sedimentBed"))
    {
        word isBed(dict_.lookup("sedimentBed"));
        if (isBed=="on")
        {
            bedExist_ = true;
        }
        else if (isBed=="off")
        {
            bedExist_ = false;
            checkMeshMotion_();
        }
        else
        {
            FatalError << "wrong keyword wrong keyword!"
                << " possible options are: on off" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    else
    {
        Info << "no sediment bed" << endl;
        bedExist_ = false;
    }
}

void Foam::sedbedManager::checkFaMeshOrientation_() const
{
    if (not bedExist_)
    {
        return;
    }
    forAll(aMesh_->faceLabels(), i)
    {
        double res = g_.value() & aMesh_->faceAreaNormals()[i];
        if (res <= 0)
        {
            FatalError
                << "sediment bed and gravity orientation "
                << "are not consistent" << endl
                << "g and patch faces dot product "
                << "should be positive" << endl;
            Info << abort(FatalError) << endl;
        }
    }
}

void Foam::sedbedManager::getPatchesID()
{
    List<word> bedPatchesNames(dict_.lookup("sedimentBedPatches"));
    forAll(bedPatchesNames, i)
    {
        word patchName = bedPatchesNames[i];
        label patchID = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchID==-1)
        {
            FatalError << "bedPatch " << patchName
                << " does not exist" << endl
                << "existing patches are:" << endl
                << mesh_.boundaryMesh().names() << endl;
            Info << abort(FatalError) << endl;
        }
        bedPatchesNames_.append(patchName);
        bedPatchesID_.append(patchID);
    }
}

void Foam::sedbedManager::checkMeshMotion_()
{
    if (dict_.found("meshMotion"))
    {
        word meshMotionState(dict_.lookup("meshMotion"));
        if (meshMotionState=="on")
        {
            meshMotion_ = true;
        }
        else if (meshMotionState=="off")
        {
            meshMotion_ = false;
            Info << "no mesh motion, "
                << "bedload and shields are computed "
                << "but bed will not move" << endl;
        }
        else
        {
            FatalError << "wrong keyword! possible options are: on off,"
                << " default is on" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    else
    {
        meshMotion_ = true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
