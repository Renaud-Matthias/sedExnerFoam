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

#include "bedloadModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModel::bedloadModel()
:
    bedloadModel("Meyer-Peter")
{}

Foam::bedloadModel::bedloadModel
(
    const word& modelName
)
:
    critShields0_(dimless, 0.047)
{
    modelName_ = modelName;
    setCoefs(modelName_);
}

Foam::bedloadModel::bedloadModel
(
    const word& modelName,
    const dimensionedScalar critShields0
)
:
    bedloadModel(modelName)
{
    critShields0_ = critShields0;
}

Foam::bedloadModel::bedloadModel
(
    scalar alpha,
    scalar aExp,
    scalar bExp,
    const dimensionedScalar& critShields0
)
:
    bedloadModel()
{
    modelName_ = "custom";
    critShields0_ = critShields0;
    setCoefs(alpha, aExp, bExp);
}

Foam::bedloadModel::bedloadModel(const dictionary& dict)
:
    bedloadModel()
{
    if (dict.found("bedloadFormula"))
    {
        word name(dict.lookup("bedloadFormula"));
        modelName_ = name;

        if (name=="custom")
        {
            scalar alpha = readScalar(dict.lookup("alpha"));
            scalar aExp = readScalar(dict.lookup("aExp"));
            scalar bExp = readScalar(dict.lookup("bExp"));
            setCoefs(alpha, aExp, bExp);
        }
        else
        {
            setCoefs(name);
        }
    }

    if (dict.found("criticalShields"))
    {
        scalar Shc = readScalar(dict.lookup("criticalShields"));
        critShields0_.value() = Shc;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModel::~bedloadModel()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::bedloadModel::setCoefs
(
    scalar alpha,
    scalar aExp,
    scalar bExp
)
{
    alpha_ = alpha;
    aExp_ = aExp;
    bExp_ = bExp;
}

void Foam::bedloadModel::setCoefs(const word& modelName)
{
    if (modelName=="Meyer-Peter")
    {
        alpha_ = 8;
        aExp_ = 0;
        bExp_ = 1.5;
    }
    else if (modelName=="Nielsen")
    {
        alpha_ = 12;
        aExp_ = 0.5;
        bExp_ = 1;
    }
    else
    {
        FatalError << "bedloadModel " << modelName
            << " for transport formula not a valid entry" << endl
            << "possible names are: Meyer-Peter,  Nielsen" << endl;
        Info << abort(FatalError) << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::bedloadModel::criticalShields0() const
{
    return critShields0_;
}

Foam::tmp<Foam::areaVectorField> Foam::bedloadModel::qb
(
    const areaVectorField& shields,
    const areaScalarField& critShields,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& g,
    const dimensionedScalar& dS
) const
{
    dimensionedScalar einsteinNumber =
        Foam::sqrt(((rhoS/rhoF) - 1) * g * Foam::pow(dS, 3));

    // in case qb is a vector, is direction is needed
    dimensionedScalar smallVal(dimless, SMALL);

    return alpha_ * einsteinNumber
        * (shields / (mag(shields) + smallVal))
        * Foam::pow(Foam::mag(shields), aExp_)
        * Foam::pow
        (
            pos(Foam::mag(shields) - critShields)
            * (Foam::mag(shields) - critShields), bExp_
        );
}

void Foam::bedloadModel::output(Ostream& os) const
{
    os << "bedload model type: " << modelName_ << endl
        << "qb = alpha sqrt((s-1)gd^3) shields^a (shields-critShields)^b"
        << endl << "alpha: " << alpha_ << endl
        << "a: " << aExp_ << endl
        << "b: " << bExp_ << endl
        << "critShields: " << critShields0_.value() << endl;
}

Foam::Ostream& Foam::operator<<(Ostream& os, const bedloadModel& model)
{
    model.output(os);
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
