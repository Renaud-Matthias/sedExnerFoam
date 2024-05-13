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

#include "sedimentBed.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedimentBed::sedimentBed
(
    dictionary& dict,
    Foam::fvMesh& mesh,
    const meshObjects::gravity& g
)
:
    bedExist_(false),
    bedMotion_(false),
    avalanche_(false),
    dict_(dict),
    mesh_(mesh),
    g_(g)
{
    checkBedExistence_();
    if (bedExist_)
    {
        aMesh_.reset(new faMesh(mesh_));
        vsm.reset(new volSurfaceMapping(aMesh_));
        getPatchesID();
        checkFaMeshOrientation_();
        
        Info << "faMesh patches ID : " << bedPatchesID_ << endl;
        Info << "faMesh patches names : " << bedPatchesNames_ << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedimentBed::~sedimentBed()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// public member functions

const faMesh& Foam::sedimentBed::aMesh()
{
    return aMesh_.ref();
}

bool Foam::sedimentBed::isAvalanche()
{
    if (avalanche_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

labelList Foam::sedimentBed::bedPatchesID()
{
    return bedPatchesID_.clone();
}

bool Foam::sedimentBed::exist() const
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

bool Foam::sedimentBed::bedMotion() const
{
    if (bedMotion_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Protected member functions

void Foam::sedimentBed::checkBedExistence_()
{
    word isBed(dict_.lookup("bedload"));
    if (isBed=="on")
    {
        bedExist_ = true;
        checkBedMotion_();
        checkAvalancheModel_();
    }
    else if (isBed=="off")
    {
        bedExist_ = false;
    }
    else
    {
        FatalError << "wrong keyword for entry bedload!"
            << " possible options are: on off" << endl;
        Info << abort(FatalError) << endl;
    }
}

void Foam::sedimentBed::checkAvalancheModel_()
{
    if (dict_.found("avalanche"))
    {
        word avModel(dict_.lookup("avalanche"));
        if (avModel=="on")
        {
            avalanche_ = true;
        }
        else if (avModel=="off")
        {
            avalanche_ = false;
        }
        else
        {
            FatalError << "wrong keyword for entry avalanche,"
                << " possible options are: on off" << endl;
            Info << abort(FatalError) << endl;
        }
    }
}

void Foam::sedimentBed::checkFaMeshOrientation_() const
{
    if (not bedExist_)
    {
        return;
    }
    Info << "g : " << g_.value() << endl;
    forAll(aMesh_->faceLabels(), facei)
    {
        double res = g_.value() & aMesh_->faceAreaNormals()[facei];
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

void Foam::sedimentBed::getPatchesID()
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

void Foam::sedimentBed::checkBedMotion_()
{
    if (dict_.found("bedMotion"))
    {
        word bedMotionState(dict_.lookup("bedMotion"));
        if (bedMotionState=="on")
        {
            bedMotion_ = true;
            Info << "bed motion activated" << endl;
        }
        else if (bedMotionState=="off")
        {
            bedMotion_ = false;
            Info << "no bed motion, "
                << "bedload and shields are computed "
                << "but bed will not move" << endl;
        }
        else
        {
            FatalError << "wrong keyword for entry bedMotion!"
                << " possible options are: on off,"
                << " default is on" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    else
    {
        bedMotion_ = true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
