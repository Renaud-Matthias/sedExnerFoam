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
#ifndef avalancheVinent_H
#define avalancheVinent_H

#include "dictionary.H"
#include "dimensionedTypes.H"
#include "areaFields.H"


namespace Foam
{
namespace bedloadModels
{

class avalancheVinent
{
protected:

    //- dictionary containing parameters fro bedloadModel
    const dictionary& dict_;
    
    //- true if avalanche is on, default is true
    bool avalanche_;

    //- Set avalanche on/off from dict_
    void setAvalanche();

    //- maximum possible bedload related to avalanche
    dimensionedScalar Qav_;

public:

    //- Constructor
    avalancheVinent
    (
        const dictionary& dict
    );

    //- Destructor
    ~avalancheVinent();

    //- Compute avalanche related bedload
    tmp<vectorField> avalanche
    (
        const scalarField& beta,
        const vectorField& slopeDir,
        const scalar& betaRep
    ) const;
    
};
}  // end namespace bedloadModels
}  // end namespace Foam

#endif
