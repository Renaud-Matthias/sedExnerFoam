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
#ifndef HindranceModel_H
#define HindranceModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
namespace settlingModels
{

class HindranceModel
{
protected:
    const dictionary& dict_;
    
    public:
    
    //- Runtime type information
    TypeName("hindranceModel");
    
    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HindranceModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
    
    // Constructor
    HindranceModel(const dictionary& dict);

    // Selectors
    static autoPtr<HindranceModel> New
    (
        const dictionary& dict
    );

    //- Destructor
    virtual ~HindranceModel();

    // Member functions
    virtual tmp<volScalarField> hindrance
    (
        const volScalarField& C
    ) const = 0;
};
} // End namespace Foam
} // End namespace settlingModels

#endif
