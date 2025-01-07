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
#include "Camenen.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bedloadModels
{
    defineTypeNameAndDebug(Camenen, 0);
    addToRunTimeSelectionTable(bedloadModel, Camenen, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModels::Camenen::Camenen(const dictionary& dict)
:
    bedloadModel(dict),
    alpha_(12),
    aExp_(1.5),
    bCoef_(4.5)
{
     Info << "bedload model type: Camenen" << endl
         << "qb* = 8 (shields-critShields)^1.5" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModels::Camenen::~Camenen()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField> Foam::bedloadModels::Camenen::qb
(
    const areaVectorField& shields,
    const areaScalarField& critShields,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& g,
    const dimensionedScalar& dS
) const
{
    // compute einstein number
    dimensionedScalar numEin =
        einsteinNumber(rhoS, rhoF, g, dS);

    dimensionedScalar smallVal(dimless, SMALL);

    return alpha_ * numEin
        * (shields / (mag(shields) + smallVal))
        * Foam::pow(Foam::mag(shields), aExp_)
        * Foam::exp(bCoef_ * Foam::mag(shields) / critShields);
}
