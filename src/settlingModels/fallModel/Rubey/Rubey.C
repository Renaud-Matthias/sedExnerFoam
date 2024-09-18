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
#include "Rubey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(Rubey, 0);
    addToRunTimeSelectionTable(fallModel, Rubey, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::Rubey::Rubey(const dictionary& dict)
:
    fallModel(dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::Rubey::~Rubey()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::settlingModels::Rubey::Ufall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar dstar = Dstar(dS, rhoS, rhoF, nuF, g);
    dimensionedScalar s = rhoS/rhoF;
    dimensionedScalar Ws = Foam::sqrt(dS*g*(s-1)) * (
        Foam::sqrt(0.667 + 36/Foam::pow(dstar, 3.))
        - Foam::sqrt(36/Foam::pow(dstar, 3.))
    );
    return Ws;
}
