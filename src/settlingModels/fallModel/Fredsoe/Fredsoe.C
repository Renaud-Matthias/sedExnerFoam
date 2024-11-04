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
#include "Fredsoe.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace settlingModels
{
    defineTypeNameAndDebug(Fredsoe, 0);
    addToRunTimeSelectionTable(fallModel, Fredsoe, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::settlingModels::Fredsoe::Fredsoe(const dictionary& dict)
:
    fallModel(dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::settlingModels::Fredsoe::~Fredsoe()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::settlingModels::Fredsoe::Ufall0
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

Foam::dimensionedScalar Foam::settlingModels::Fredsoe::nextValue
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

Foam::dimensionedScalar Foam::settlingModels::Fredsoe::solveUfall0
(
    const dimensionedScalar& dS,
    const dimensionedScalar& s,
    const dimensionedScalar& nuF,
    const dimensionedScalar& g
) const
{
    dimensionedScalar Ws(dimVelocity, 0.1);
    dimensionedScalar err(dimless, 1);
    dimensionedScalar eps(dimless, 1e-3);
    dimensionedScalar Wsnext(nextValue(Ws, dS, s, nuF, g));
    label i(0);
    label imax(1000);  // maximum number of iteration
    while (i<imax and err>eps)
        {
            i += 1;
            Ws = Wsnext;
            Wsnext = nextValue(Ws, dS, s, nuF, g);
            err = Foam::mag((Ws - Wsnext)/Wsnext);
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
