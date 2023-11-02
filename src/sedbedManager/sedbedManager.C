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
    bedExist_(false),
    meshMotion_(false),
    avalanche_(false),
    dict_(dict),
    mesh_(mesh),
    g_(g),
    diffFilter_(dimensionSet(0, 2, 0, 0, 0, 0, 0), 0)
{
    checkBedExistence_();
    if (bedExist_)
    {
        aMesh_.reset(new faMesh(mesh_));
        vsm.reset(new volSurfaceMapping(aMesh()));
        getPatchesID();
        checkFaMeshOrientation_();
        checkExnerFilter_();
        
        Info << "faMesh patches ID : " << bedPatchesID_ << endl;
        Info << "faMesh patches names : " << bedPatchesNames_ << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedbedManager::~sedbedManager()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// public member functions

faMesh& Foam::sedbedManager::aMesh()
{
    return aMesh_.ref();
}

areaVectorField& Foam::sedbedManager::qb()
{
    return qb_.ref();
}

areaVectorField& Foam::sedbedManager::shields()
{
    return shields_.ref();
}

areaVectorField& Foam::sedbedManager::zb()
{
    return zb_.ref();
}

areaVectorField& Foam::sedbedManager::dH()
{
    return dH_.ref();
}

edgeScalarField& Foam::sedbedManager::phiqb()
{
    return phiqb_.ref();
}

areaScalarField& Foam::sedbedManager::angleBedSlope()
{
    return beta_.ref();
}

bool Foam::sedbedManager::isAvalanche()
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

void Foam::sedbedManager::filterExner()
{
    if (filterExner_)
    {
        // inverse distance bewteen faces centers
        edgeScalarField deltaCoeff =
            aMesh().edgeInterpolation::deltaCoeffs();

        // project zb along vertical direction
        areaScalarField zbProj = zb() & (-g_/ mag(g_));
        zb() = zb() - zbProj * (-g_/ mag(g_));
        // compute diffusivity for laplacian operator
        diffFilter_ = alphaFilter_ * max(1/pow(deltaCoeff, 2));
        
        // predictor step
        zbPred.ref() = zbProj
            + fac::laplacian(diffFilter_, zbProj);

        // corrector step
        for (label icorr=0; icorr<Ncorr_; icorr++)
        {
            dzCorr.ref() = zbProj - zbPred.ref();
            dzCorr.ref() = dzCorr.ref()
                + fac::laplacian(diffFilter_, dzCorr.ref());
            zbProj = zbPred.ref() + dzCorr.ref();
        }
        zb() = zb() + zbProj * (-g_/ mag(g_));
        zb().correctBoundaryConditions();
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
            checkMeshMotion_();
            checkAvalancheModel_();
        }
        else if (isBed=="off")
        {
            bedExist_ = false;
        }
        else
        {
            FatalError << "wrong keyword"
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

void Foam::sedbedManager::checkAvalancheModel_()
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

void Foam::sedbedManager::checkExnerFilter_()
{
    word switchFilter(dict_.lookup("filterExner"));
    if (switchFilter=="on")
    {
        filterExner_ = true;
        alphaFilter_ = readScalar(dict_.lookup("alphaFilter"));
        Ncorr_ = readLabel(dict_.lookup("Ncorrection"));
    }
    else if (switchFilter=="off")
    {
        filterExner_ = false;
    }
    else
    {
        FatalError << "wrong keyword for entry filterExner,"
            << " possible options are: on off" << endl;
        Info << abort(FatalError) << endl;
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
            FatalError << "wrong keyword for entry meshMotion!"
                << " possible options are: on off,"
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
