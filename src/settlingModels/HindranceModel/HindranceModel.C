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
#include "HindranceModel.H"


namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(HindranceModel, 0);
    defineRunTimeSelectionTable(HindranceModel, dictionary);
}
}

// Constructor

Foam::settlingModels::HindranceModel::HindranceModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// Destructor

Foam::settlingModels::HindranceModel::~HindranceModel()
{}

// Selector

Foam::autoPtr<Foam::settlingModels::HindranceModel>
Foam::settlingModels::HindranceModel::New
(
    const dictionary& dict
)
{
    word hindranceModelType(dict.get<word>("type"));

    Info<< "Selecting hindranceModel "
        << hindranceModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(hindranceModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "HindranceModel::New(const dictionary&) : " << endl
            << "    unknown hindranceModelType type "
            << hindranceModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid hindranceModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<HindranceModel>(cstrIter()(dict));
}
