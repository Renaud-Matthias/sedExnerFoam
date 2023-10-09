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


edgeScalarField deltaCoeff =
    aMesh.edgeInterpolation::deltaCoeffs();

// diffusivity for laplacian operator
diffFilter = alphaFilter * min(1/pow(deltaCoeff, 2));
Info << "filter diffusivity : " << diffFilter.value() << endl;
// predictor step
Zbpred = Zb + fac::laplacian(diffFilter, Zb);

for (label icorr=0; icorr<Ncorr; icorr++)
{
    dzCorr = Zb - Zbpred;
    dzCorr = dzCorr + fac::laplacian(diffFilter, dzCorr);
    Zb = Zbpred + dzCorr;
}

Zb.correctBoundaryConditions();
