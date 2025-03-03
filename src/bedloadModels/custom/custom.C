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
#include "custom.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bedloadModels
{
    defineTypeNameAndDebug(custom, 0);
    addToRunTimeSelectionTable(bedloadModel, custom, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModels::custom::custom(const dictionary& dict)
:
    bedloadModel(dict),
    alpha_(dict_.get<scalar>("alpha")),
    aExp_(dict_.get<scalar>("aExp")),
    bExp_(dict_.get<scalar>("bExp"))
{
     Info << "bedload model type: custom" << endl
         << "qb* = " << alpha_
         << " shields^" << aExp_
         << " (shields-critShields)^" << bExp_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModels::custom::~custom()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField> Foam::bedloadModels::custom::qb
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

    dimensionedScalar smallVal(dimless, SMALL);

    return alpha_ * numEin
        * (shields / (mag(shields) + smallVal))
        * Foam::pow(Foam::mag(shields), aExp_)
        * Foam::pow
        (
            pos(Foam::mag(shields) - critShields)
            * (Foam::mag(shields) - critShields),
            bExp_
        );
}
