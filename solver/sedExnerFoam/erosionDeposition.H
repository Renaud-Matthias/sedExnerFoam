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

// Access to horizontaly projected variables
faMeshProjection& meshProj = aMeshProjPtr.ref();  // projected area mesh
// area of faces projected on horizontal plane
const scalarField& SfProj = meshProj.S();

// reference to erosion and deposition terms
volScalarField& Ero = erosionPtr.ref();
volScalarField& Dep = depositionPtr.ref();

volScalarField& Cs = CsPtr.ref();  // suspended load
volVectorField& Ws = WsPtr.ref();  // settling velocity

areaScalarField& DE = DEptr.ref();  // deposition - erosion
surfaceScalarField& bedPhi = bedPhiPtr.ref();  // mesh flux through faces

const faMesh& aMesh = bed.aMesh();

for (const label patchID : bed.bedPatchesID())
{
    // deposition boundary field
    scalarField& DepPatch = Dep.boundaryFieldRef()[patchID];
    // erosion boundary field
    scalarField& EroPatch = Ero.boundaryFieldRef()[patchID];
    // reference to bed fvPatch
    const fvPatch& bedPatch = mesh.boundary()[patchID];
    const scalarField& bedPhiPatch = bedPhi.boundaryFieldRef()[patchID];
    
    //const scalarField& diffSedPatch = diffSed.boundaryFieldRef()[patchID];
    // get gradient of Cs on bed patch faces

    forAll(bedPatch, facei)
    {
        // label of owner cell of face
        label idNCell = bedPatch.faceCells()[facei];
        scalar ci = Cs[idNCell];
        vector wsi = Ws[idNCell];
        scalar Sface = aMesh.S()[facei];
        scalar Vcell = mesh.V()[idNCell];
        scalar bedPhiFace = bedPhiPatch[facei];
        const vector& nFace = aMesh.faceAreaNormals()[facei];
        // deposition due to settling
        DepPatch[facei] = ci * ((wsi & nFace) - bedPhiFace/Sface);
        // deposition due to diffusion
        //DepPatch[facei] +=
        DE[facei] = DepPatch[facei] - EroPatch[facei] * (Vcell/Sface);
    }
}

/*DE.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface<scalar>
    (
        Dep.boundaryFieldRef()
    )
  - bed.vsm.ref().mapToSurface<scalar>
    (
        Ero.boundaryFieldRef()
    );*/

// deposition and erosion contribution in Exner equation
forAll(bed.aMesh().areaCentres(), facei)
{
    dVsed[facei] += DE[facei] * SfProj[facei];
}
