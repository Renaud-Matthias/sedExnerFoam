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

#include "settlingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModel::settlingModel
(
    const dictionary& dict
)
:
    dict_(dict)
{
        Info << "Initialization of settlingModel" << endl;
        fallDict_ = dict_.subDict("fallModel");
        Info << "Initialization of fallModel" << endl;
        fallModel_ = settlingModels::fallModel::New (fallDict_);
        hindranceDict_ = dict_.subDict("hindranceModel");
        Info << "Initialization of hindranceModel" << endl;
        hindranceModel_ = settlingModels::hindranceModel::New (hindranceDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModel::~settlingModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::settlingModel::Ufall
(
    const volScalarField& C,
    const dimensionedScalar& Cmax,
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar ufall0(fallModel_->getUfall0(dS, rhoS, rhoF, nuF, g));
    return ufall0*hindranceModel_->hindrance(C, Cmax);
}
