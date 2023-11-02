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

word patchName = "bed";

label aPatchID = mesh.boundaryMesh().findPatchID(patchName);

//Info << "faMesh faces label: " << aMesh.faceLabels() << endl;
//Info << "boundaryPatch faces labell: " << mesh.boundaryMesh()[aPatchID].

// set-up interpolator
primitivePatchInterpolation patchInterpolator(mesh.boundaryMesh()[aPatchID]);

// compute face displacements
forAll(aMesh.areaCentres(), facei)
{
    scalar dh = Zb[facei] - aMesh.areaCentres()[facei].z();
    deltaH[facei].z() = dh;
    // show faces a and b new coordinates
}

// interpolate dh on mesh points

vectorField dispVals = patchInterpolator.faceToPointInterpolate(deltaH);

/*Info << ", min(|deltaH|) : " << min(mag(dispVals));
Info << ", max(|deltaH|) : " << max(mag(dispVals)) << endl;
Info << ", mean(deltaH) : " << average(deltaH.internalFieldRef()) << endl;
*/

pointField newPoints = mesh.points();

const labelList& bedPointsLabels = mesh.boundaryMesh()[patchName].meshPoints();
const pointField& lp = mesh.boundaryMesh()[patchName].localPoints();

forAll(bedPointsLabels, i)
{
    label pointPatch = bedPointsLabels[i];
    //newPoints[pointPatch] += dispVals[i];
    newPoints[pointPatch] = lp[i] + dispVals[i];
}

Info << "mesh movePoints" << endl;
mesh.movePoints(newPoints);
