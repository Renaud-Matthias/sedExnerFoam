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

#include "avalancheVinent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedloadModels::avalancheVinent::avalancheVinent
(
    const dictionary& dict
)
:
    dict_(dict),
    Qav_(dimVelocity * dimLength, 0.01)
{
    setAvalanche();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bedloadModels::avalancheVinent::~avalancheVinent()
{}

// * * * * * * * * * * * * * * * Member funcions * * * * * * * * * * * * * * //

void Foam::bedloadModels::avalancheVinent::setAvalanche()
{
    Qav_.value() = dict_.lookupOrDefault<scalar>("Qav", 1e-3);
    if (Qav_.value() < 0)
    {
        FatalError
            << "Qav value must be positive" << endl;
        Info << abort(FatalError) << endl;
    }
}

Foam::tmp<Foam::vectorField>
Foam::bedloadModels::avalancheVinent::avalanche
(
    const scalarField& beta,
    const vectorField& slopeDir,
    const scalar& betaRep
) const
{
    return Qav_.value() * slopeDir * Foam::pos(beta - betaRep)
        * (
            Foam::tanh(Foam::tan(beta))
            - Foam::tanh(Foam::tan(betaRep))
        )
        / (1 - Foam::tanh(Foam::tan(betaRep)));
}
