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
#include "Miedema.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalShieldsModels
{
    defineTypeNameAndDebug(Miedema, 0);
    addToRunTimeSelectionTable(criticalShieldsModel, Miedema, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Miedema::Miedema(const dictionary& dict)
:
    criticalShieldsModel(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Miedema::~Miedema()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::criticalShieldsModels::Miedema::criticalShields0
(
    const dimensionedScalar& Dstar
) const
{
    dimensionedScalar critShields =
        (0.2285 / Foam::pow(Dstar, 1.02))
        + 0.0575 * (1 - Foam::exp(-0.0225*Dstar));
    return critShields;
}
