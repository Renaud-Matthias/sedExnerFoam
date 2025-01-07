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


areaVectorField& shields = shieldsPtr.ref();

// compute wallShearStress on each patch belonging to sediment bed
for (const label patchID : bed.bedPatchesID())
{
    // compute wall shear stress on patch
    vectorField& ssp = shieldsVf.boundaryFieldRef()[patchID];
    const vectorField& Sfp = mesh.Sf().boundaryField()[patchID];
    const scalarField& magSfp =
        mesh.magSf().boundaryField()[patchID];
    volSymmTensorField Reff =
        rhoF * turbulence->nuEff() * dev(twoSymm(fvc::grad(U)));
    const symmTensorField& Reffp = Reff.boundaryField()[patchID];
    
    // compute shields number
    ssp = ((-Sfp / magSfp) & Reffp)
        / (mag(g).value() * dS.value()
        * (rhoS.value() - rhoF.value()));
    
    // remove bed normal component, shear stress only
    ssp += (ssp & (-Sfp / magSfp)) * (Sfp / magSfp);
}

// map volumic shields number to area shields number
shields.primitiveFieldRef() =
    bed.vsm.ref().mapToSurface<vector>
    (
        shieldsVf.boundaryFieldRef()
    );

// remove bed normal component from shields, shear stress only
//shields -=
//    (shields & bed.aMesh().faceAreaNormals())
//    * (shields / (Foam::mag(shields) + SMALL));
