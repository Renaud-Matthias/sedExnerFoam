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

#include "criticalShieldsModel.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalShieldsModels
{
    defineTypeNameAndDebug(criticalShieldsModel, 0);
    defineRunTimeSelectionTable(criticalShieldsModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::criticalShieldsModel::criticalShieldsModel
(
    const dictionary& dict
)
:
    dict_(dict),
    slopeCorrection_(false)
    //critShields0_(dimensionedScalar(dimless, 0.047))
{
    if (dict_.found("slopeCorrection"))
    {
        word switchSlopeCorr = word(dict_.lookup("slopeCorrection"));
        if (switchSlopeCorr == "on")
        {
            slopeCorrection_ = true;
        }
        else if (switchSlopeCorr != "off")
        {
            FatalError << "wrong value for entry slopeCorrection, "
                << "possible options are: on off" << endl;
        Info << abort(FatalError) << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::criticalShieldsModels::criticalShieldsModel::~criticalShieldsModel()
{}

// * * * * * * * * * * * * * * * * Selectors *  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::criticalShieldsModels::criticalShieldsModel>
Foam::criticalShieldsModels::criticalShieldsModel::New
(
    const dictionary& dict
)
{
    word criticalShieldsModelType(dict.get<word>("type"));

    Info<< "Selecting critical Shields number Model "
        << criticalShieldsModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(criticalShieldsModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "criticalShieldsModel::New(const dictionary&) : " << endl
            << "    unknown criticalShieldsModelType type "
            << criticalShieldsModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid criticalShieldsModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<criticalShieldsModel>(cstrIter()(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::criticalShieldsModels::criticalShieldsModel::slopeCorrection
(
    areaScalarField& critShields,
    const areaVectorField& shields,
    const scalarField& slopeAngle,
    const vectorField& slopeDir,
    scalar reposeAngle
) const
{
    if (slopeCorrection_)
    {
        // direction of shear stress on bed
        vectorField stressDir = shields / (Foam::mag(shields) + SMALL);
        forAll(critShields, facei)
        {
            scalar beta = slopeAngle[facei];
            scalar alpha = Foam::acos(slopeDir[facei] & stressDir[facei]);
            scalar betap = Foam::min(beta, reposeAngle);
            critShields[facei] *=
                (
                    Foam::cos(betap)
                    * Foam::sqrt
                    (
                        1 - Foam::pow
                        (
                            Foam::sin(alpha)
                            * Foam::tan(betap)
                            * (1/Foam::tan(reposeAngle)),
                            2
                        )
                    )
                    - Foam::cos(alpha)
                    * Foam::sin(betap) / Foam::tan(reposeAngle)
                );
        }
    }
}

Foam::dimensionedScalar
Foam::criticalShieldsModels::criticalShieldsModel::viscousDiameter
(
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& nu,
    const dimensionedScalar& g
)
{
    dimensionedScalar s(rhoS/rhoF);  // density ratio
    // compute viscous diameter
    dimensionedScalar dv = Foam::pow(Foam::pow(nu, 2)/(g*(s-1)), 1/3);
    return dS / dv;
}
