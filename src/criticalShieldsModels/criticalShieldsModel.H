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
#ifndef criticalShieldsModel_H
#define criticalShieldsModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
namespace criticalShieldsModels
{

class criticalShieldsModel
{

    protected:

    //- dictionary containing parameters fro criticalShieldsModel
    const dictionary& dict_;


    public:
    
    //- Runtime type information
    TypeName("criticalShieldsModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        criticalShieldsModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    // Constructor
    criticalShieldsModel
    (
        const dictionary& dict
    );

    // Selectors
    static autoPtr<criticalShieldsModel> New
    (
        const dictionary& dict
    );

    //- Destructor
    virtual ~criticalShieldsModel();


    // Member functions

    //- Compute viscous diameter
    dimensionedScalar viscousDiameter
    (
        const dimensionedScalar& dS,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nu,
        const dimensionedScalar& g
    );

    //- Return critical Shields number, computed from dimless diameter
    virtual dimensionedScalar criticalShields0
    (
        const dimensionedScalar& Dstar
    ) const = 0;
    
};
}  // End namespace Foam
}  // End namespace criticalShieldsModels

#endif
