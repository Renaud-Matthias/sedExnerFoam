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

const faMesh& aMesh = bed.aMesh();

areaScalarField& dH = dHptr.ref();
areaScalarField& dVsed = dVsedPtr.ref();

// faMesh horizontal projection
const projectedFaMesh& meshProj = bed.aProjMesh();
// area of projected mesh faces
const scalarField& SfProj = meshProj.S();

#include "balanceBedload.H"

Info << "bedload balanced" << endl;

//#include "bedFormVelocity.H"

Info << "bedform velocity computed" << endl;

if (switchSuspension=="on")
{
    #include "erosionDeposition.H"
}

scalar deltaT = runTime.time().deltaT().value();  // current time step

// solve exner equation explicitly with Adams-Bashforth scheme
if (ABorder==0)
{
    // Adams-Bashforth order 0, euler explicit scheme
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoAccFactor
            * (1/CsMax.value())
            * (deltaT / SfProj[facei])
            * dVsed[facei];
    }
}
else if (ABorder==1)
{
    scalar deltaT0 = runTime.time().deltaT0().value();  // old time step

    scalar coefT0 = deltaT * (1 + 0.5 * (deltaT / deltaT0));
    scalar coefT00 = - 0.5 * deltaT * deltaT / deltaT0;

    // Adams-Bashforth order 1
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoAccFactor
            * (1/CsMax.value())
            * (coefT0 * dVsed[facei] + coefT00 * dVsed.oldTime()[facei])
            * (1 / SfProj[facei]);
    }
}

if (filterExner)
{
    scalar dl = 0;
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        scalar newdl = Foam::mag(
            xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]]);
        if (newdl>dl)
        {
            dl = newdl;
        }
    }
    
    scalar Kex = Foam::pow(dl, 2) * 0.1;
    label NfiltZb = 3;

    for (label j=0; j<NfiltZb; j++)
    {
        scalarField lapldH(bed.aMesh().nFaces());

        for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
        {
            vector Xon =
                xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]];
            //scalar w = (ne & Xen) / (ne & Xon);
            scalar graddHe =
                (
                    (
                        dH[neighbour[edgei]]
                        - dH[owner[edgei]]
                    ) * (Xon & LeProj[edgei])
                ) / Foam::pow(Foam::mag(Xon), 2);

            lapldH[owner[edgei]] =
                (Kex / SfProj[owner[edgei]]) * graddHe;
            lapldH[neighbour[edgei]] =
                -(Kex / SfProj[neighbour[edgei]]) * graddHe;
        }
        dH.primitiveFieldRef() += lapldH;
    }
}

if (bed.bedMotion())
{
    #include "meshMove.H"
}
