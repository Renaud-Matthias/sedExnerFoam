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
#ifndef FallModel_H
#define FallModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
namespace settlingModels
{

class fallModel
{
protected:

    //- Dictionary containing parameters for fallModel
    const dictionary& dict_;

    //- Pointer falling velocity of a lone particle
    mutable dimensionedScalar* ufallPtr_;

    // Protected member functions
    
    //- Falling velocity of a lone particle in a fluid
    virtual dimensionedScalar Ufall0
    (
        const dimensionedScalar& dS,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nuF,
        const dimensionedScalar& g
    ) const = 0;
    
    public:
    
    //- Runtime type information
    TypeName("fallModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fallModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
    
    // Constructors

    //- Construct from dictionary
    fallModel(const dictionary& dict);

    // Selectors
    
    static autoPtr<fallModel> New
    (
        const dictionary& dict
    );

    // Destructors
    
    virtual ~fallModel();

    // Member functions

    //- Get falling velocity of a lone particle in a fluid
    dimensionedScalar getUfall0
    (
        const dimensionedScalar& dS,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nuF,
        const dimensionedScalar& g
    ) const;
};

} // End namespace settlingModels
} // End namespace Foam

#endif
