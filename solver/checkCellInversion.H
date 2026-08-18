/*---------------------------------------------------------------------------* \
Copyright (C) 2026 Arthur Lopez, Matthias Renaud, Cyrille Bonamy,
                   Julien Chauchat and contributors

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

const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

// TODO, allow different bed patch names
const label bedPatchID = bMesh.findPatchID("bed");

const polyPatch& bedPatch = bMesh[bedPatchID];

boolList isBedPoint(mesh.nPoints(), false);

const labelList& bedMeshPoints = bedPatch.meshPoints();
forAll(bedMeshPoints, i)
{
    isBedPoint[bedMeshPoints[i]] = true;
}

const labelList& bedFaceCells = bedPatch.faceCells();

vector upDir = -eg;

const pointField& pts = mesh.points();
const labelListList& cellPointAddr = mesh.cellPoints();

label nInverted = 0;
scalar maxPenetration = 0.0;

forAll(bedFaceCells, i)
{
    const label cellI = bedFaceCells[i];
    const labelList& cellNodes = cellPointAddr[cellI];

    scalar maxBedHeight = -GREAT;
    scalar minOtherHeight = GREAT;

    forAll(cellNodes, j)
    {
        const label pointI = cellNodes[j];
        const scalar h = pts[pointI] & upDir;

        if (isBedPoint[pointI])
        {
            maxBedHeight = max(maxBedHeight, h);
        }
        else
        {
            minOtherHeight = min(minOtherHeight, h);
        }
    }

    if (maxBedHeight > -GREAT && minOtherHeight < GREAT)
    {
        if (maxBedHeight > minOtherHeight)
        {
            nInverted++;
            const scalar penetration = maxBedHeight - minOtherHeight;
            maxPenetration = max(maxPenetration, penetration);

        }
    }
}

reduce(nInverted, sumOp<label>());
reduce(maxPenetration, maxOp<scalar>());

if (nInverted > 0)
{
    WarningInFunction
        << "checkCellInversion.H: detected " << nInverted
        << " inverted cell(s) on patch bed (max penetration = "
        << maxPenetration << " m)."
        << endl;
}
