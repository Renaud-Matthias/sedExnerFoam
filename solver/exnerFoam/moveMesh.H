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

#include "filterExner.H"

//Info << aMesh.areaCentres().internalField() << endl;
//Info << Zb << endl;

word patchName = "bed";

label aPatchID = mesh.boundaryMesh().findPatchID(patchName);

// set-up interpolator
primitivePatchInterpolation patchInterpolator(mesh.boundaryMesh()[aPatchID]);

// compute face displacements
forAll(aMesh.areaCentres(), facei)
{
    scalar dh = Zb[facei] - aMesh.areaCentres()[facei].z();
    deltaH[facei].z() = dh;
}

// interpolate dh on mesh points

vectorField dispVals = patchInterpolator.faceToPointInterpolate(deltaH);

Info << ", min(|deltaH|) : " << min(mag(dispVals));
Info << ", max(|deltaH|) : " << max(mag(dispVals)) << endl;
Info << ", mean(deltaH) : " << average(deltaH.internalFieldRef()) << endl;

/*pointField newPoints = mesh.points();

labelList bedPointsLabels = mesh.boundaryMesh()[patchName].meshPoints();

forAll(bedPointsLabels, i)
{
    label pointPatch = bedPointsLabels[i];
    newPoints[pointPatch] += dispVals[i];
}

mesh.movePoints(newPoints);*/

pointDisplacement.boundaryFieldRef()[aPatchID] == dispVals;

bool meshChanged = mesh.update();

if (meshChanged)
{
    Info << "mesh has been updated" << endl;
}
