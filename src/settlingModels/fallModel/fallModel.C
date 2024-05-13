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
#include "fallModel.H"


namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(fallModel, 0);
    defineRunTimeSelectionTable(fallModel, dictionary);
}
}

// Constructor

Foam::settlingModels::fallModel::fallModel
(
    const dictionary& dict
)
:
    dict_(dict),
    ufallPtr_(nullptr)
{}


// Destructor

Foam::settlingModels::fallModel::~fallModel()
{}

// Selector

Foam::autoPtr<Foam::settlingModels::fallModel>
Foam::settlingModels::fallModel::New
(
    const dictionary& dict
)
{
    word fallModelType(dict.get<word>("type"));

    Info<< "Selecting fallModel "
        << fallModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(fallModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "fallModel::New(const dictionary&) : " << endl
            << "    unknown fallModelType type "
            << fallModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid fallModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<fallModel>(cstrIter()(dict));
}

// Member functions

Foam::dimensionedScalar Foam::settlingModels::fallModel::getUfall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    if (ufallPtr_==nullptr)
    {
        ufallPtr_ = new dimensionedScalar(Ufall0(dS, rhoS, rhoF, nuF, g));
        
        dimensionedScalar& ufall = *ufallPtr_;
        Info << "falling velocity of a lone particle: "
            << ufall.value() << " m/s" << endl;
    }

    return *ufallPtr_;
}
