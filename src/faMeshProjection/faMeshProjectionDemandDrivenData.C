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

#include "faMeshProjection.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// vector a is projected on a plane normal to vector b
vector projectNormalPlane
(
    const vector& a,
    const vector& planeNormal
)
{
    vector aProj = a - (a & planeNormal) * (planeNormal / mag(planeNormal));
    return aProj;
}

}  // End of namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshProjection::calcLe() const
{
    if (LePtr_)
    {
        FatalError
            << "LePtr_ already allocated" << endl;
        Info << abort(FatalError);
    }

    LePtr_ = new vectorField(nEdges_);

    const vectorField& Le = mesh_.Le();

    vectorField& LeProj = *LePtr_;

    forAll(Le, edgei)
    {
        LeProj[edgei] = projectNormalPlane(Le[edgei], projectNormal_);
    }
}

void Foam::faMeshProjection::calcAreaCentres() const
{
    if (areaCentresPtr_)
    {
        FatalError
            << "faceCentresPtr_ already allocated" << endl;
        Info << abort(FatalError);
    }
    
    areaCentresPtr_ = new vectorField(nFaces_);
    // face centers coordinates on curved faMesh
    const vectorField& areaCentres =
        mesh_.areaCentres().internalField();

    vectorField& areaCentresProj = *areaCentresPtr_;

    forAll(areaCentres, facei)
    {
        areaCentresProj[facei] =
            projectNormalPlane(areaCentres[facei], projectNormal_);
    }
}

void Foam::faMeshProjection::calcEdgeCentres() const
{
    if (edgeCentresPtr_)
    {
        FatalError
            << "edgeCentresPtr_ already allocated" << endl;
        Info << abort(FatalError);
    }

    edgeCentresPtr_ = new vectorField(nEdges_);
    // face centers coordinates on curved faMesh
    const vectorField& edgeCentres =
        mesh_.edgeCentres().internalField();

    vectorField& edgeCentresProj = *edgeCentresPtr_;

    forAll(edgeCentres, edgei)
    {
        edgeCentresProj[edgei] =
            projectNormalPlane(edgeCentres[edgei], projectNormal_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
