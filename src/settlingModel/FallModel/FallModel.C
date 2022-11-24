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
#include "FallModel.H"


namespace Foam
{
    defineTypeNameAndDebug(FallModel, 0);
    defineRunTimeSelectionTable(FallModel, dictionary);
}

// Constructor

Foam::FallModel::FallModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


// Destructor

Foam::FallModel::~FallModel()
{}

// Selector

Foam::autoPtr<Foam::FallModel> Foam::FallModel::New
(
    const dictionary& dict
)
{
    word fallModelType(dict.get<word>("fallModel"));

    Info<< "Selecting fallModel "
        << fallModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(fallModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "FallModel::New(const dictionary&) : " << endl
            << "    unknown fallModelType type "
            << fallModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid fallModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<FallModel>(cstrIter()(dict));
}
