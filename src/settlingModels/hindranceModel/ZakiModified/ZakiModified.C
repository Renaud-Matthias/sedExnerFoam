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
#include "ZakiModified.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(ZakiModified, 0);
    addToRunTimeSelectionTable(hindranceModel, ZakiModified, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::ZakiModified::ZakiModified(const dictionary& dict)
:
    hindranceModel(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::ZakiModified::~ZakiModified()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::settlingModels::ZakiModified::hindrance
(
    const volScalarField& C,
    const dimensionedScalar& Cmax
) const
{
    scalar nExponent(dict_.get<scalar>("n"));
    return pow(1 - C, nExponent - 1)
        * pow((1 - C / Cmax) * neg(C - Cmax), Cmax);
}
