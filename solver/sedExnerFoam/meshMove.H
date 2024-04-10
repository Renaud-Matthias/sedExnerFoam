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

// projected area mesh
faMeshProjection& meshProj = aMeshProjPtr.ref();

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

vectorField& dispVals = dispValsPtr.ref();
scalarField& dHpoints = dHpointsPtr.ref();
dHpoints = Zero;

// get vertices displacement through interpolation for each bed patch
for (const label patchID : bed.bedPatchesID())
{
    // Access to face centers in projected horizontal plane
    const vectorField& xFacesProj = meshProj.areaCentres();

    // Access to points coordinates in projected horizontal plane
    const vectorField& xPointsProj = meshProj.pointCoords();

    // normalise interpolation weight
    scalarField normWeight(aMesh.nPoints());
    normWeight = Zero;

    // gravity oriented unit vector
    vector eg = (g / Foam::mag(g)).value();

    // loop over all faces
    for (label facei=0; facei < aMesh.nFaces(); facei++)
    {
        const face& f = aMesh.faces()[facei];
        const vector& xf = xFacesProj[facei];
        // number of vertices defining face f
        const label nVerts = f.size();
        
        // loop over points belonging to face
        for (label pi=0; pi < nVerts; pi++)
        {
            // label of point pi
            const label piGlob = f.thisLabel(pi);

            scalar areaPF = 0;

            // coordinates of vertices
            const vector& xv = xPointsProj[piGlob];
            const vector& xn = xPointsProj[f.nextLabel(pi)];
            const vector& xp = xPointsProj[f.prevLabel(pi)];

            // vector from vertice to face center
            vector vp = xf - xv;
            // vector from vertice to right edge center
            vector vredg = 0.5 * (xn - xv);
            // vector from vertice to left edge center
            vector vledg = 0.5 * (xp - xv);

            areaPF +=
                0.5 * mag(vredg ^ vp)
                + 0.5 * mag(vp ^ vledg);

            dHpoints[piGlob] += areaPF * dH[facei];

            normWeight[piGlob] += areaPF;
        }
    }
    // Normalise weight of interpolation
    dHpoints /= normWeight;

    dispVals -= eg * dHpoints;

    // interpolate displacement on points, vertical orientation
    //dispVals -= patchInterp.faceToPointInterpolate(dH) * eg;
    
    Info << "mesh vertices motion, max:" << max(mag(dispVals));
    Info << " m , min:" << min(mag(dispVals));
    Info << " m , mean:" << average(mag(dispVals)) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
}

// reset some fields to 0
areaScalarField& dVsed = dVsedPtr.ref();
forAll(bed.aMesh().areaCentres(), facei)
{
    dVsed[facei] = 0;
}
