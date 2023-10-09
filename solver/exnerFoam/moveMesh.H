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

// global label
label pointabGlob(168);
label faceaGlob(1668);
label facebGlob(1669);

// faMesh local labels for faces a and b
label faceaLoc(-1);
label facebLoc(-1);

// local point ab
label pointabLoc(-1);

forAll(aMesh.faceLabels(), i)
{
    if (aMesh.faceLabels()[i]==faceaGlob)
    {
        faceaLoc = i;
    }
    else if (aMesh.faceLabels()[i]==facebGlob)
    {
        facebLoc = i;
    }
}

Info << "global labels" <<endl;
Info << "face a: " << faceaGlob << " " << mesh.Cf()[faceaGlob] << endl;
Info << "point ab: " << pointabGlob << " " << mesh.points()[pointabGlob] << endl;
Info << "face b: " << facebGlob << " " << mesh.Cf()[facebGlob] << endl;

Info << "local faMesh faces labels" << endl;
Info << "face a: " << faceaLoc << " " << aMesh.areaCentres()[faceaLoc] << endl;
Info << "face b: " << facebLoc << " " << aMesh.areaCentres()[facebLoc] << endl;

// set-up interpolator
primitivePatchInterpolation patchInterpolator(mesh.boundaryMesh()[aPatchID]);

// compute face displacements
forAll(aMesh.areaCentres(), facei)
{
    scalar dh = Zb[facei] - aMesh.areaCentres()[facei].z();
    deltaH[facei].z() = dh;
    // show faces a and b new coordinates
    if (facei==faceaLoc)
    {
        Info << "dh face a: " << aMesh.areaCentres()[facei] << " dh: " << deltaH[facei] << endl;
    }
    else if (facei==facebLoc)
    {
        Info << "dh face b: " << aMesh.areaCentres()[facei] << " dh: " << deltaH[facei] << endl;
    }
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
    if (pointPatch==pointabGlob)
    {
        pointabLoc = i;
    }
}

Info << "dh point ab: " << dispVals[pointabLoc] << endl;
Info << "expected new point ab: " << newPoints[pointabGlob] << endl;

Info << "mesh movePoints" << endl;
mesh.movePoints(newPoints);

Info << "point ab after mesh update: " << mesh.points()[pointabGlob] << endl;
