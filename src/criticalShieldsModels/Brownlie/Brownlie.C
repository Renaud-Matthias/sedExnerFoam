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
#include "Brownlie.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalShieldsModels
{
    defineTypeNameAndDebug(Brownlie, 0);
    addToRunTimeSelectionTable(criticalShieldsModel, Brownlie, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Brownlie::Brownlie(const dictionary& dict)
:
    criticalShieldsModel(dict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::Brownlie::~Brownlie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::criticalShieldsModels::Brownlie::calcCriticalShields0
(
    const dimensionedScalar& Dstar
) const
{
    if (critShields0_)
    {
        FatalError
            << "critShields0_ already allocated" << endl;
    }
    critShields0_.reset
        (
            new dimensionedScalar
            (
                "critShields0",
                dimless,
                Zero
            )
        );
    dimensionedScalar& critShields0 = critShields0_.ref();
    critShields0 =
        (0.22 / Foam::pow(Dstar, 0.9))
        + Foam::pow(
            0.0610,
            -7.7 * Foam::pow(Dstar, -0.9)
        );
}
