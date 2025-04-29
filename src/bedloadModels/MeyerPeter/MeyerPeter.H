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
#ifndef MeyerPeter_H
#define MeyerPeter_H

#include "bedloadModel.H"

namespace Foam
{
namespace bedloadModels
{

class MeyerPeter
:
    public bedloadModel
{
protected:

    scalar alpha_;

    scalar bExp_;

public:

    // Runtime type information
    TypeName("MeyerPeter");

    //- Construct from components
    MeyerPeter(const dictionary& dict);


    //- Destructor
    virtual ~MeyerPeter();


    // Member functions

    //- Compute bedload from Shields and sediment properties
    tmp<areaVectorField> getQb
    (
        const areaVectorField& shields,
        const areaScalarField& critShields,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& g,
        const dimensionedScalar& dS,
        const dimensionedScalar& nuF
    ) const;
};

} // End namespace bedloadModels
} // End namespace Foam

#endif
