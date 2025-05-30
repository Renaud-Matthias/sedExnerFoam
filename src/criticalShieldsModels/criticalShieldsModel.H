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
#include "areaFields.H"
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

    //- pointer to base critical Shields number
    mutable autoPtr<dimensionedScalar> critShields0_;

    //- pointer to dimless diameter
    autoPtr<dimensionedScalar> Dstar_;

    //- true if slope correction activated
    bool slopeCorrection_;

    //- Compute base critical SHields number
    virtual void calcCriticalShields0
    (
        const dimensionedScalar& Dstar
    ) const = 0;


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

    //- Apply slope correction on critical Shields
    void slopeCorrection
    (
        areaScalarField& critShields,
        const vectorField& shields,
        const scalarField& slopeAngle,
        const vectorField& slopeDir,
        scalar reposeAngle
    ) const;

    //- Compute viscous diameter
    dimensionedScalar viscousDiameter
    (
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nu,
        const dimensionedScalar& g
    );

    //- Compute dimless diameter
    dimensionedScalar Dstar
    (
        const dimensionedScalar& dS,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nu,
        const dimensionedScalar& g
    );

    //- Return critical Shields number, computed from dimless diameter
    const dimensionedScalar& criticalShields0
    (
        const dimensionedScalar& Dstar
    );
};
}  // End namespace Foam
}  // End namespace criticalShieldsModels

#endif
