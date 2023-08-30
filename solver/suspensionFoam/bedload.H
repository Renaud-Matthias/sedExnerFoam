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

// compute wallShearStress on each patch belonging to sediment bed
for (const label patchID : faBed.bedPatchesID())
{
    Info << "patch ID : " << patchID << endl;
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
}

// map volumic shields number to area shields number
faBed.shields.ref().primitiveFieldRef() =
    faBed.vsm.ref().mapToSurface<vector>
    (
        shieldsVf.boundaryFieldRef()
    );
// compute bedload with Meyer Peter Muller law
// no threshold at the moment
dimensionedScalar mpmE(dimless, 8);

Info << "compute qb using MPM law" << endl;

dimensionedScalar smallVal(dimless, SMALL);
faBed.qb.ref() =
    mpmE * (faBed.shields.ref() / (smallVal + mag(faBed.shields.ref())))
    * Foam::sqrt((((rhoS / rhoF) - 1) * mag(g) * Foam::pow(dS, 3)))
    * Foam::pow(pos(mag(faBed.shields.ref()) - critShields)
    * (mag(faBed.shields.ref()) - critShields), 1.5);


// compute qb flux through fa mesh edges
Info << "combute flux of qb through faces edges" << endl;
faBed.phib.ref() =
    linearEdgeInterpolate(faBed.qb.ref()) & faBed.aMesh_.ref().Le();


// solve exner explicitly (euler explicit)
Info << "explicit exner solve" << endl;
faBed.zb.ref() = faBed.zb.ref() + runTime.time().deltaT()
    * vector(0, 0, 1) * ( -fac::edgeIntegrate(faBed.phib.ref()));

// compute faces displacements
Info << "compute faces displacements" << endl;
forAll(faBed.aMesh_.ref().areaCentres(), facei)
{
    vector dh =
        faBed.zb.ref()[facei] - faBed.aMesh_.ref().areaCentres()[facei];
    faBed.deltaH.ref()[facei] = dh;
}

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

if (faBed.meshMotion())
{
        #include "meshMove.H"
}
