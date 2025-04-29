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
#include "vanRijn.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bedloadModels
{
    defineTypeNameAndDebug(vanRijn, 0);
    addToRunTimeSelectionTable(bedloadModel, vanRijn, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModels::vanRijn::vanRijn(const dictionary& dict)
:
    bedloadModel(dict),
    alpha_(0.053),
    aExp_(2.1),
    bExp_(-0.3)
{
     Info << "bedload model type: vanRijn" << endl
         << "qb* = 0.053 (shields/critShields - 1)**2.1 / Dstar**0.3" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModels::vanRijn::~vanRijn()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField> Foam::bedloadModels::vanRijn::getQb
(
    const areaVectorField& shields,
    const areaScalarField& critShields,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& g,
    const dimensionedScalar& dS,
    const dimensionedScalar& nuF
) const
{
    // compute einstein number
    dimensionedScalar numEin =
        einsteinNumber(rhoS, rhoF, g, dS);

    // compute dimless diameter
    dimensionedScalar Dstar = dS * Foam::pow(
        ((rhoS/rhoF - 1) * g) / Foam::pow(nuF, 2),
        1.0/3.0
    );

    dimensionedScalar smallVal(dimless, SMALL);

    return alpha_ * numEin
        * (shields / (mag(shields) + smallVal))
        * Foam::pow
        (
            (Foam::mag(shields)/(critShields + smallVal) - 1)
            * Foam::pos(Foam::mag(shields)/(critShields + smallVal) - 1),
            aExp_
        )
        * Foam::pow(Dstar, bExp_);
}
