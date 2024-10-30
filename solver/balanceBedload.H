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


// Access to bedload variable
areaVectorField& qb = qbPtr.ref();

// project qb on horizontal plane
vectorField qbProj = meshProj.project(qb);

// projected edge centers coordinates
const vectorField& xEdgesProj = meshProj.edgeCentres();
// projected face centers coordinates
const vectorField& xFacesProj = meshProj.areaCentres();
// projected edge vectors, oriented toward owner face
const vectorField& LeProj = meshProj.Le();
// face owner list
const labelList& owner = meshProj.edgeOwner();
// face neighbour list
const labelList& neighbour = meshProj.edgeNeighbour();

// get qb on edges, linear interpolation
for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
{
    // distance between edgei center and owner face center
    scalar distOwner =
        Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        // distance between edgei center and neighbour center
    scalar distNeighbour =
        Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
    // linear interpolate qb on edge i
    vector qbEdgei =
        (
            distOwner * qbProj[owner[edgei]]
            + distNeighbour * qbProj[neighbour[edgei]]
        )
        * (1 / (distOwner + distNeighbour));

    scalar phiqbEdgei = qbEdgei & LeProj[edgei];
    
    dVsed[owner[edgei]] -= phiqbEdgei;
    
    dVsed[neighbour[edgei]] += phiqbEdgei;
}

// get qb value for edges on boundary patches
forAll(aMesh.boundary(), patchi)
{
    const faPatchField<Vector<double>>& QbBound =
        qb.boundaryField()[patchi];

    if (aMesh.boundary()[patchi].type()!="empty")
    {
        // faMesh edge adressing
        label edgei = aMesh.boundary()[patchi].start();

        forAll(aMesh.boundary()[patchi], patchEdgei)
        {
            const vector edgeLeProj =
                meshProj.project
                (
                    aMesh.Le().boundaryField()[patchi][patchEdgei]
                );

            vector qbEdgei = meshProj.project(QbBound[patchEdgei]);

            scalar phiqbEdgei = qbEdgei & edgeLeProj;

            // apply boundary condition to exner equation
            dVsed[owner[edgei]] -= phiqbEdgei;

            edgei++;
        }
    }
}
