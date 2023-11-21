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

    scalar critShields_;
    
public:

    //const static List<word> modelNames = {"Nielsen"};

    // Constructors

    //- Construct from name of model
    bedloadModel
    (
        const word& modelName
    );

    //- Construct from name of model and critical shields number
    bedloadModel
    (
        const word& modelName,
        scalar critShields
    );

    //- Construct from coefficients
    bedloadModel
    (
        scalar alpha,
        scalar a,
        scalar b,
        scalar critShields
    );

    // Destructor
    virtual ~bedloadModel();

    // Member functions

    //- Return bedload value associated to shields values
    template<class Type>
    Type qb
    (
        const Type& shields,
        const scalar& rhoS,
        const scalar& rhoF,
        const scalar& g,
        const scalar& dS
    )
    {
        scalar einsteinNumber =
            Foam::sqrt(((rhoS/rhoF) - 1) * g * Foam::pow(dS, 3));

        Type relShields = shields - critShields_;

        return alpha_ * einsteinNumber
            * Foam::pow(shields, aExp_)
            * Foam::pow(pos(relShields) * (relShields), bExp_);
    }

    //- Write bedloadModel to Ostream
    void output(Ostream& os) const;
    
};

Ostream& operator<<(Ostream& os, const bedloadModel& model);

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
