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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bedloadModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModel::bedloadModel
(
    const word& modelName
)
{
    if (modelName=="Meyer-Peter")
    {
        alpha_ = 8;
        aExp_ = 0;
        bExp_ = 1.5;;
    }
    else if (modelName=="Nielsen")
    {
        alpha_ = 12;
        aExp_ = 0.5;
        bExp_ = 1;
    }
    else
    {
        Info << "not implemented yet" << endl;
    }

    // default critical Shields number value
    critShields_ = 0.047;

    Info << "model " << modelName
        << ", alpha: " << alpha_
            << ", a: " << aExp_
            << ", b: " << bExp_
            << ", critShields: " << critShields_ << endl;
}

Foam::bedloadModel::bedloadModel
(
    const word& modelName,
    scalar critShields
)
: bedloadModel(modelName)
{
    critShields_ = critShields;
}

Foam::bedloadModel::bedloadModel
(
    scalar alpha,
    scalar a,
    scalar b,
    scalar critShields
)
:
    alpha_(alpha),
    aExp_(a),
    bExp_(b),
    critShields_(critShields)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModel::~bedloadModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bedloadModel::output(Ostream& os) const
{
    os << "bedload model type: " << modelName_ << endl
        << "qb = alpha sqrt((s-1)gd^3) shields^a (shields-critShields)^b"
        << endl << "alpha: " << alpha_
        << ", a: " << aExp_
        << ", b: " << bExp_
        << ", critShields: " << critShields_ << endl;
}

Foam::Ostream& Foam::operator<<(Ostream& os, const bedloadModel& model)
{
    model.output(os);
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
