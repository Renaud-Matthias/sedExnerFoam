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
#ifndef ShieldsFixedValue_H
#define ShieldsFixedValue_H

#include "criticalShieldsModel.H"

namespace Foam
{
namespace criticalShieldsModels
{

class fixedValue
:
    public criticalShieldsModel
{

    public:

    // Runtime type information
    TypeName("fixedValue");


    // Constructors

    //- Construct from components
    fixedValue(const dictionary& dict);


    // Destructor
    virtual ~fixedValue();


    // Member functions

    //- Return critical Shields number, computed from dimless diameter
    dimensionedScalar criticalShields0
    (
        const dimensionedScalar& Dstar
    ) const;
};

} // End namespace criticalShieldsModels
} // End namespace Foam

#endif
