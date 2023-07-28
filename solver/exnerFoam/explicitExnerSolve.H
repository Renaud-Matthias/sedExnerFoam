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

// update bedload flux through edges
phib = linearEdgeInterpolate(Qb) & aMesh.Le();

rhsExnerEq = fac::edgeIntegrate(phib);

// euler explicit
if (iterTimeLoop == 1)
{
    Zb = Zb + runTime.time().deltaT() * (- rhsExnerEq);
}

// Adam Bashforth 1 (order 2)
else
{
    Zb = Zb + 0.5 * runTime.time().deltaT() *
        (
            3 * rhsExnerEq
            - rhsExnerEq.oldTime()
        );
}

/*Zb += runTime.time().deltaT() *
    (
        - fac::edgeIntegrate(phib)
        //fac::div(Qb)  // bedload
        //- fac::laplacian(Da, Zb)  // avalanche related diffusion
    );*/
