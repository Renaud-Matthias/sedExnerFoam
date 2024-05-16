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
#include "Soulsby.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalShieldsModels
{
    defineTypeNameAndDebug(Soulsby, 0);
    addToRunTimeSelectionTable(criticalShieldsModel, Soulsby, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Soulsby::Soulsby(const dictionary& dict)
:
    criticalShieldsModel(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Soulsby::~Soulsby()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::criticalShieldsModels::Soulsby::criticalShields0
(
    const dimensionedScalar& Dstar
) const
{
    dimensionedScalar critShields =
        (0.3 / (1 + 1.2*Dstar))
        + 0.055 * (1 - Foam::exp(0.02*Dstar));
    return critShields;
}
