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

// predictor step
forAll(aMesh.areaCentres(), facei)
{
    if (facei!=0 and facei!=aMesh.nFaces()-1)
    {
        Zbpred[facei] = (1-2*alphaFilter) * Zb[facei]
            + alphaFilter * (Zb[facei-1] + Zb[facei+1]);
    }
    else
    {
        Zbpred[facei] = Zb[facei];
    }
}
for (label icorr=0; icorr<Ncorr; icorr++)
{
    scalarField dzbCorr = Zb - Zbpred;
    forAll(aMesh.areaCentres(), facei)
    {
        if (facei!=0 and facei!=aMesh.nFaces()-1)
        {
            dzbCorr = (1-2*alphaFilter) * dzbCorr[facei]
                + alphaFilter * (dzbCorr[facei-1] + dzbCorr[facei+1]);
            Zb[facei] = Zbpred[facei] + dzbCorr[facei];
        }
    }
}

Zb.correctBoundaryConditions();
