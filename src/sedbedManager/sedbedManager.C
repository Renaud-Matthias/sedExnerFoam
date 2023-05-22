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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    findBedPatches();
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

void Foam::sedbedManager::findBedPatches()
{
    if (dict_.found("sedimentBedPatches"))
    {
        List<word> bedPatchesNames(dict_.lookup("sedimentBedPatches"));
        forAll(bedPatchesNames, i)
        {
            word patchName = bedPatchesNames[i];
            label patchID = checkPatchExistence(patchName);
            const polyPatch patch = mesh_.boundaryMesh()[patchID];
            checkPatchOrientation(patch);
            bedPatchesNames_.append(patchName);
            bedPatchesID_.append(patchID);
        }
        bedExist_ = true;
    }
    else
    {
        Info << "no sediment bed in the domain" << endl;
        bedExist_ = false;
    }
}

Foam::List<Foam::word> Foam::sedbedManager::getPatchesNames() const
{
    return bedPatchesNames_.clone();
}

Foam::label Foam::sedbedManager::checkPatchExistence(word patchName) const
{
    label patchID = mesh_.boundaryMesh().findPatchID(patchName);
    if (patchID==-1)
    {
        FatalError
            << "bedPatch " << patchName
            << "does not exist" << endl
            << "existing patches are:" << endl
            << mesh_.boundaryMesh().names() << endl;
        Info << abort(FatalError) << endl;
    }
    return patchID;
}

void Foam::sedbedManager::checkPatchOrientation
(
    const polyPatch& patch
) const
{
    forAll(patch.faceAreas(), i)
    {
        double res = g_.value() & patch.faceAreas()[i];
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
