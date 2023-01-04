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

#ifndef SettlingModel_H
#define SettlingModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "HindranceModel.H"
#include "FallModel.H"

namespace Foam
{

class SettlingModel
{
protected:

    const dictionary& dict_;

    dictionary fallDict_;

    dictionary hindranceDict_;

    autoPtr<settlingModels::FallModel> FallModel_;

    autoPtr<settlingModels::HindranceModel> HindranceModel_;
    
public:

    //- Constructor
    SettlingModel
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~SettlingModel();

    //- Member functions
    tmp<Foam::volScalarField> Ufall
    (
        const volScalarField& C,
        const dimensionedScalar& dS,
        const dimensionedScalar& rhoS,
        const dimensionedScalar& rhoF,
        const dimensionedScalar& nuF,
        const dimensionedScalar& g
    ) const;
};

} // End namespace foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
