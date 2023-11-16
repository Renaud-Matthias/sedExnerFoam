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


// depositon source term
volScalarField& C = CPtr.ref();  // suspended load
surfaceScalarField& phip = phipPtr.ref();  // flux(U+Ws)
areaScalarField& depC = depositionPtr.ref();
// rate of deposition in m/s
depC.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface(phip)
    * bed.vsm.ref().mapToSurface(C.boundaryField());

depC.primitiveFieldRef() =
    depC.primitiveFieldRef()
    / (bed.aMesh().S() * (bed.aMesh().faceAreaNormals() & (gUnit)));

// solve exner explicitly
// deposition contribution
dH.primitiveFieldRef() +=
    runTime.time().deltaT() * (1/Cmax) * depC;
