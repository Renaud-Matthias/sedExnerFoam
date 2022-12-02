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
#include "Fixed.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(Fixed, 0);
    addToRunTimeSelectionTable(FallModel, Fixed, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::Fixed::Fixed(const dictionary& dict)
:
    FallModel(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::Fixed::~Fixed()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::settlingModels::Fixed::Ufall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS
) const
{
    Info << dict_ << endl;
    Info << "check 2" << endl;
    //scalar UfallValue(dict_.get<scalar>("value"));
    Info << dict_ << endl;
    word UfallType(dict_.get<word>("type"));
    scalar UfallValue(0.3);
    Info << UfallType << endl;
    Info << "check 3" << endl;
    return UfallValue;
}
