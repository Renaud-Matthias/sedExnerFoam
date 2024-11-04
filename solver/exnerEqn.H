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
faMeshProjection& meshProj = aMeshProjPtr.ref();  // projected area mesh
// area of projected mesh faces
const scalarField& SfProj = meshProj.S();

if (switchBedload=="on")
{
    #include "balanceBedload.H"
}

if (switchSuspension=="on")
{
    #include "erosionDeposition.H"
}

scalar deltaT = runTime.time().deltaT().value();

// solve exner equation explicitly (euler explicit)
forAll(aMesh.areaCentres(), facei)
{
    dH[facei] =
        (deltaT / SfProj[facei])
        * (1/CsMax.value())
        * dVsed[facei];
}

if (bed.bedMotion())
{
    #include "meshMove.H"
}
