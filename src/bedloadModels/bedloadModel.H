/*---------------------------------------------------------------------------* \
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

#ifndef bedloadModel_H
#define bedloadModel_H

#include "dictionary.H"
#include "areaFields.H"
#include "dimensionedTypes.H"


namespace Foam
{

class bedloadModel
{
protected:

    // bedload model
    // qbs = alpha * shields^aExp - (shields - shieldsC)^bExp

    word modelName_;

    scalar alpha_;

    scalar aExp_;

    scalar bExp_;

    dimensionedScalar critShields0_;

    void setCoefs(scalar alpha, scalar aExp, scalar bExp);

    void setCoefs(const word& modelName);
    
public:

    // Constructors

    //- Default constructor
    bedloadModel();

    //- Construct from name of model
    bedloadModel
    (
        const word& modelName
    );

    //- Construct from name of model and base critical Shields number
    bedloadModel
    (
        const word& modelName,
        dimensionedScalar critShields0
    );

    //- Construct from coefficients
    bedloadModel
    (
        scalar alpha,
        scalar aExp,
        scalar bExp,
        const dimensionedScalar& critShields0
    );

    //- Construct from dictionary
    bedloadModel
    (
        const dictionary& dict
    );


    // Destructor
    
    virtual ~bedloadModel();


    // Member functions

    //- Return constant reference to base critical Shields number
    const dimensionedScalar& criticalShields0() const;

    //- Compute bedload from base Shields and sediment properties
    template<class Type>
    Type qb
    (
        const Type& shields,
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
        Type qbDirection =
            shields / (mag(shields) + smallVal);

        Type qbOut = alpha_ * einsteinNumber
            * qbDirection
            * Foam::pow(Foam::mag(shields), aExp_)
            * Foam::pow
            (
                pos(Foam::mag(shields) - critShields0_)
                * (Foam::mag(shields) - critShields0_), bExp_
            );

        return qbOut;
    }

    //- Compute bedload from Shields and sediment properties
    tmp<areaVectorField> qb
    (
        const areaVectorField& shields,
        const areaScalarField& critShields,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& g,
        const dimensionedScalar& dS
    ) const;

    //- Write bedloadModel to Ostream
    void output(Ostream& os) const;
    
};

Ostream& operator<<(Ostream& os, const bedloadModel& model);

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
