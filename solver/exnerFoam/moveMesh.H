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

//Info << aMesh.areaCentres().internalField() << endl;
//Info << Zb << endl;

word patchName = "bed";

label aPatchID = mesh.boundaryMesh().findPatchID(patchName);

primitivePatchInterpolation patchInterpolator(mesh.boundaryMesh()[aPatchID]);

// predictor step
Zbedges = linearEdgeInterpolate(Zb);
//patchInterpolator.faceToEdgeInterpolate(Zb);
//scalarField Zbpred(Zb.clone());

// filter zb
// compute face displacements
diffFilter.primitiveFieldRef() = alphaFilter * aMesh.S();

for (label iout=0; iout<Nouter; iout++)
{
    Zbpred = Zb + fac::laplacian(diffFilter, Zb);
    dzCorr = Zb - Zbpred;
    for (label iin=0; iin<Ninner; iin++)
    {
        dzCorr = dzCorr + fac::laplacian(diffFilter, dzCorr);
    }
    Zb = Zbpred + dzCorr;
    Zb.correctBoundaryConditions();
}

/*
for (label iout=0; iout<Nouter; iout++)
{
    // predictor step
    forAll(aMesh.areaCentres(), i)
    {
        if (i!=0 and i!=aMesh.nFaces()-1)
        {
            Zbpred[i] = 0.5 * Zbpred[i] + 0.25 * (Zbpred[i-1] + Zbpred[i+1]);
        }
    }
    scalarField dzbCorr = Zb - Zbpred;
    for (label iin=0; iin<Ninner; iin++)
    {
        forAll(aMesh.areaCentres(), i)
        {
            if (i!=0 and i!=aMesh.nFaces()-1)
            {
                dzbCorr = 0.5 * dzbCorr[i]
                    + 0.25 * (dzbCorr[i-1] + dzbCorr[i+1]);
            }
        }
    }
    Zbpred += dzbCorr;
}

    Zb.primitiveFieldRef() = Zbpred;*/


// compute face displacements
forAll(aMesh.areaCentres(), facei)
{
    scalar dh = Zb[facei] - aMesh.areaCentres()[facei].z();
    deltaH[facei] = dh;
}

// interpolate dh on mesh points

scalarField pointValues = patchInterpolator.faceToPointInterpolate(deltaH);

Info << ", min(deltaH) : " << min(pointValues);
Info << ", max(deltaH) : " << max(pointValues) << endl;

//Info << "pointValues : " << pointValues << endl;

// Mesh motion
//vectorField newPoints(mesh.points().size(), 0);

//Info << "delta H : " << pointValues << endl;

//Info << "boundaryMesh : " << mesh.boundaryMesh()[patchName] << endl;

//Info << "bed patch ID : " << aPatchID
//<< ", patchName : " << patchName << endl;

//Info << "boundaryMesh points coordinates"
//<< mesh.boundaryMesh()[patchName].meshPoints() << endl;

pointField newPoints = mesh.points();

labelList bedPointsLabels = mesh.boundaryMesh()[patchName].meshPoints();

forAll(bedPointsLabels, i)
{
    label pointPatch = bedPointsLabels[i];
    newPoints[pointPatch].z() += pointValues[i];
}

mesh.movePoints(newPoints);

//bool meshChanged = mesh.update();

//Info << "mesh update : " << meshChanged << endl;
