/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    suspensionFoam

Group
    grpBasicSolvers

Description
    Scalar transport and incompressible turbulent flow solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{T} + \div \left(\vec{U} T\right) - \div \left(D_T \grad T \right)
        = S_{T}
    \f]

    Where:
    \vartable
        C          | Volumic fraction
        \vec{U}    | Velocity
        \vec{R}    | Stress tensor
        p          | Pressure
        \vec{S}_U  | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        C       | Passive scalar
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "CMULES.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "SettlingModel.H"
#include "sedbedManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Scalar transport and incompressible turbulent flow solver."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFaFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();
    
    if (not LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "CEqn.H"
            
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            
            #include "UEqn.H"
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        if (faBed.exist())
        {
            //const volSymmTensorField& Reff = turbulence->devReff();
            for (const label patchID : faBed.bedPatchesID())
            {
                Info << "patch ID : " << patchID << endl;
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
                // map volumic shields number to area shields number
                faBed.shields.ref().primitiveFieldRef() =
                    faBed.vsm.ref().mapToSurface<vector>
                    (
                        shieldsVf.boundaryFieldRef()
                    );
                
                // Compute bedload with Meyer Peter Muller law
                // no threshold at the moment
                dimensionedScalar mpmE(dimless, 8);
                Info << "compute qb using MPM law" << endl;
                dimensionedScalar smallVal(dimless, SMALL);
                faBed.qb.ref() = mpmE * (faBed.shields.ref() / (smallVal + mag(faBed.shields.ref())))
                    * Foam::sqrt((((rhoS / rhoF) - 1) * mag(g) * Foam::pow(dS, 3)))
                    * Foam::pow(pos(mag(faBed.shields.ref()) - critShields)
                    * (mag(faBed.shields.ref()) - critShields), 1.5);
            }
            //Info << "Info : " << turbulence->devReff() << endl;
        }

        if (runTime.writeTime())
        {
            if (faBed.exist())
            {
                // map area<>Fields to vol<>Fields
                faBed.vsm.ref().mapToVolume
                    (
                        faBed.shields.ref(),
                        shieldsVf.boundaryFieldRef()
                    );
                faBed.vsm.ref().mapToVolume
                    (
                        faBed.qb.ref(),
                        qbVf.boundaryFieldRef()
                    );
            }
            runTime.write();

            runTime.printExecutionTime(Info);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
