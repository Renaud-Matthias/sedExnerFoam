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
#ifndef bedloadModel_H
#define bedloadModel_H

#include "dictionary.H"
#include "volFields.H"
#include "areaFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "avalancheVinent.H"

namespace Foam
{
namespace bedloadModels
{

class bedloadModel
{

protected:

    //- dictionary containing parameters fro bedloadModel
    const dictionary& dict_;

    //- model name
    word modelName_;

    //- Coefficient of proportionality to amplify/reduce Shields number
    scalar coefShields_;

    //- true if avalanche is on
    bool avalanche_;

    //- Pointer to model for avalanche
    // bedload due slope exceeding angle of repose
    autoPtr<avalancheVinent> avalancheModel_;

    //- true if slope correction on Shields is applied
    bool slopeCorr_;

    //- coefficient for slope correction
    scalar gamma_;


public:

    //- Runtime type information
    TypeName("bedloadModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        bedloadModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Constructor
    bedloadModel
    (
        const dictionary& dict
    );

    // Selectors
    static autoPtr<bedloadModel> New
    (
        const dictionary& dict
    );

    //- Destructor
    virtual ~bedloadModel();


    // Member functions

    //- Compute bedload from base Shields and sediment properties
    tmp<areaVectorField> qb
    (
        const areaVectorField& shields,
        const areaScalarField& critShields,
        const areaVectorField& gradZb,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& g,
        const dimensionedScalar& dS,
        const dimensionedScalar& nuF
    ) const;

    //- Compute bedload from base Shields and sediment properties
    virtual tmp<areaVectorField> getQb
    (
        const areaVectorField& shields,
        const areaScalarField& critShields,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& g,
        const dimensionedScalar& dS,
        const dimensionedScalar& nuF
    ) const = 0;

    // compute and return the Einstein number
    dimensionedScalar einsteinNumber
    (
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& g,
        const dimensionedScalar& dS
    ) const;

    //- Coefficient to amplify Shields number
    scalar coefShields() const;

    //- Return true if avalanche is on
    bool avalanche() const;

    //- Return true if slope correction on Shields is applied
    bool slopeCorrection() const;

    //- Return avalanche related bedload
    tmp<vectorField> qbAvalanche
    (
        const scalarField& beta,
        const vectorField& slopeDir,
        const scalar& betaRep
    ) const;
};

}  // End namespace bedloadModels

}  // End namespace Foam

#endif
