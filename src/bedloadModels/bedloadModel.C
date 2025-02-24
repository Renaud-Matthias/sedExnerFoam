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

#include "bedloadModel.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace bedloadModels
{
    defineTypeNameAndDebug(bedloadModel, 0);
    defineRunTimeSelectionTable(bedloadModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModels::bedloadModel::bedloadModel
(
    const dictionary& dict
)
:
    dict_(dict),
    coefShields_(1),
    avalanche_(true),
    avalancheModel_(nullptr)
{
    if (dict_.found("coefShields"))
    {
        coefShields_ = readScalar(dict_.lookup("coefShields"));
        Info << "Shields number amplification factor: "
            << coefShields_ << endl;
    }
    if (coefShields_ <= 0)
    {
        FatalError
            << "Shields number amplification factor"
            << " coefShields must be positive" << endl;
        Info << abort(FatalError) << endl;
    }

    // initialize or not avalancheModel
    word switchAv(dict_.lookupOrDefault<word>("avalanche", "on"));
    if (switchAv == "on")
    {
        avalanche_ = true;
    }
    else if (switchAv == "off")
    {
        avalanche_ = false;
    }
    else
    {
        FatalError
            << "wrong entry for keyword avalanche"
                << "possible options are: on off" << endl;
        Info << abort(FatalError) << endl;
    }
    if (avalanche_)
    {
        avalancheModel_.reset
            (
                new avalancheVinent(dict_)
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModels::bedloadModel::~bedloadModel()
{}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::bedloadModels::bedloadModel>
Foam::bedloadModels::bedloadModel::New
(
    const dictionary& dict
)
{
    word bedloadModelType(dict.get<word>("type"));

    Info<< "Selecting bedloadModel: "
        << bedloadModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(bedloadModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "bedloadModel::New(const dictionary&) : " << endl
            << "    unknown bedloadModelType type "
            << bedloadModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid bedloadModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<bedloadModel>(cstrIter()(dict));
}

// * * * * * * * * * * * * * * * Member funcions * * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::bedloadModels::bedloadModel::einsteinNumber
(
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& g,
    const dimensionedScalar& dS
) const
{
    return Foam::sqrt(((rhoS/rhoF) - 1) * g * Foam::pow(dS, 3));
}

Foam::scalar Foam::bedloadModels::bedloadModel::coefShields() const
{
    return coefShields_;
}

bool Foam::bedloadModels::bedloadModel::avalanche() const
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

Foam::tmp<Foam::vectorField>
Foam::bedloadModels::bedloadModel::qbAvalanche
(
    const scalarField& beta,
    const vectorField& slopeDir,
    const scalar& betaRep
) const
{
    if (avalancheModel_==nullptr)
    {
        FatalError
            << "avalanche model not initialized, "
                << "avalanche related bedload cannot be computed" << endl;
        Info << abort(FatalError) << endl;
    }
    return avalancheModel_->avalanche(beta, slopeDir, betaRep);
}
