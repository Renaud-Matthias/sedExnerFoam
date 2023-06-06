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
    bedExist_(false), dict_(dict), mesh_(mesh), g_(g)
{
    checkBedExistence();
    if (bedExist_)
    {
        aMesh_.reset(new faMesh(mesh_));
        vsm.reset(new volSurfaceMapping(aMesh_.ref()));
        getPatchesID();
        checkFaMeshOrientation();
        
        Info << "faMesh patches ID : " << bedPatchesID_ << endl;
        Info << "faMesh patches names : " << bedPatchesNames_ << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedbedManager::~sedbedManager()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


void Foam::sedbedManager::checkBedExistence()
{
    if (dict_.found("sedimentBed"))
    {
        word isBed(dict_.lookup("sedimentBed"));
        if (isBed=="yes")
        {
            bedExist_ = true;
        }
        else if (isBed=="no")
        {
            bedExist_ = false;
        }
        else
        {
            FatalError
                << "wrong keyword" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    else
    {
        Info << "no sediment bed" << endl;
        bedExist_ = false;
    }
}

labelList Foam::sedbedManager::bedPatchesID()
{
    return bedPatchesID_.clone();
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
            FatalError
                << "bedPatch " << patchName
                << " does not exist" << endl
                << "existing patches are:" << endl
                << mesh_.boundaryMesh().names() << endl;
            Info << abort(FatalError) << endl;
        }
        bedPatchesNames_.append(patchName);
        bedPatchesID_.append(patchID);
        //bedPatches_.append(patch);
    }
}

void Foam::sedbedManager::checkFaMeshOrientation() const
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
