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

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

vectorField& dispVals = dispValsPtr.ref();

// get vertices displacement through interpolation for each bed patch
for (const label patchID : bed.bedPatchesID())
{
    primitivePatchInterpolation patchInterp(mesh.boundaryMesh()[patchID]);

    // interpolate displacement on points, vertical orientation
    dispVals -= patchInterp.faceToPointInterpolate(dH) * gUnit;
    
    Info << "mesh vertices motion, max:" << max(mag(dispVals));
    Info << " , min:" << min(mag(dispVals));
    Info << " , mean:" << average(mag(dispVals)) << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
}
