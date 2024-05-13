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
#include "Deigaard.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(Deigaard, 0);
    addToRunTimeSelectionTable(fallModel, Deigaard, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::Deigaard::Deigaard(const dictionary& dict)
:
    fallModel(dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::Deigaard::~Deigaard()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::settlingModels::Deigaard::Ufall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& rhoS,
    const dimensionedScalar& rhoF,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar s(rhoS/rhoF);
    dimensionedScalar Ws
    (
        solveUfall0(dS, s, nuF, g)
    );
    return Ws;
}

Foam::dimensionedScalar Foam::settlingModels::Deigaard::nextValue
(
    const dimensionedScalar& WsOld,
    const dimensionedScalar& dS,
    const dimensionedScalar& s,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar dragCoef(1.4+(36*nuF)/(WsOld*dS));
    dimensionedScalar WsNew(sqrt((4*(s-1)*g*dS)/(3*dragCoef)));
    return WsNew;
}

Foam::dimensionedScalar Foam::settlingModels::Deigaard::solveUfall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& s,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar Ws("Ws", dimVelocity, 0.1);
    dimensionedScalar eps("eps", dimVelocity, 1e-5);
    dimensionedScalar Wsnext(nextValue(Ws, dS, s, nuF, g));
    label i(0);
    label imax(10);
    while (i<imax and mag(Ws-Wsnext)>eps)
        {
            i += 1;
            Ws = Wsnext;
            Wsnext = nextValue(Ws, dS, s, nuF, g);
        }
    Info << "solve settling velocity, ";
    if (i>=imax)
    {
        Info << "convergence not reached : ";
        Info << i << " iterations" << endl;
    }
    else
    {
        Info << "convergence reached : ";
        Info << i << " iterations" << endl;
    }
    return Wsnext;
}
