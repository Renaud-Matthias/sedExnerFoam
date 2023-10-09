/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    exnerFoam

Description
    finiteArea equation solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt(Zb) + \div(Qb) = \laplacian(Da, Zb)
    \f]

    Where:
    \vartable
        Zb      | Bed elevationPassive scalar
        phib    | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        Zb      | bed elevation [m]
        U       | Velocity [m/s]
        Qb      | Bedload flux [m2/s]
    \endplaintable

Author
    Matthias Renaud, LEGI

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "meshTools.H"
#include "faCFD.H"
#include "motionSolver.H"
#include "primitivePatchInterpolation.H"
#include "unitConversion.H"
#include "pointMesh.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Exner equation solver."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFaFields.H"

    // set initial condition for dune tutorial
    // initiate Zb level, dune or triangle for slide

    #include "setInitialCondition.H"
    
    #include "createVolFields.H"

    Info << "\nStarting time loop\n" << endl;

    label iterTimeLoop = 0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.value() << endl;

        iterTimeLoop += 1;

        if (isMeshMoving)
        {
            theta.primitiveFieldRef() =
                Foam::acos(aMesh.faceAreaNormals() & (g / mag(g)).value());
        }
        else
        {
            theta = Foam::atan(Foam::mag(fac::grad(Zb)));
        }

        // finite area normals
        const vectorField& faNormals = aMesh.faceAreaNormals();
        
        if (bedload=="on")
        {
            forAll(aMesh.areaCentres(), i)
            {
                scalar qb = alpha *
                    Foam::pow(
                        Qwater.value().x() / (Hwater.value() - Zb[i]), beta);
                Qb[i] = vector(qb, 0, 0);
            }
        }
        
        if (avalanche=="on")
        {
            #include "updateDa.H"
        }

        Info << "max(Da) : " << max(Da) << endl;

        // bed slope correction
        forAll(faNormals, i)
        {
            vector surfaceNormal = faNormals[i];
            scalar slopeCorr = surfaceNormal & (g / mag(g)).value();
            Qb[i] /= slopeCorr;
            Da[i] /= slopeCorr;
        }
        Qb.correctBoundaryConditions();

        dimensionedScalar smallLength(dimLength, SMALL);
        Ub = Qb / (Zb + smallLength);
        
        // flux of qb through edges
        phiqb = linearEdgeInterpolate(Qb) & aMesh.Le();
        
        if (eqRes=="custom")
        {
                #include "explicitExnerSolve.H"
        }
        // with matrix construction
        else if (eqRes=="explicit")
        {
            Info << "exner equation explicit solve" << endl;
            faScalarMatrix ZbEqn
                (
                    fam::ddt(Zb)
                    + fac::div(Qb)  // explicit
                    - fam::laplacian(Da, Zb)
                );
            ZbEqn.solve();
        }
        else if (eqRes=="implicit")
        {
            phiub = linearEdgeInterpolate(Ub) & aMesh.Le();
            Info << "exner equation implicit solve" << endl;
            faScalarMatrix ZbEqn
                (
                    fam::ddt(Zb)
                    + fam::div(phiub, Zb) // implicit
                    - fac::laplacian(Da, Zb)
                );
            ZbEqn.solve();
        }

        Zb.correctBoundaryConditions();
        Qb.correctBoundaryConditions();

        if (filterType=="laplacian")
        {
            #include "filterDiff.H"
        }
        else if (filterType=="direct")
        {
            #include "filterDirect.H"
        }
        
        Info<< "bed elevation = "
            << Zb.weightedAverage(aMesh.S()).value()
            << "  Min(Zb) = " << gMin(Zb)
            << "  Max(Zb) = " << gMax(Zb)
            << endl;

        if (isMeshMoving)
        {
            #include "moveMesh.H"
        }

        if (runTime.writeTime())
        {
            vsm.mapToVolume(Zb, Zbvf.boundaryFieldRef());
            vsm.mapToVolume(Qb, Qbvf.boundaryFieldRef());
            vsm.mapToVolume(deltaH, deltaHvf.boundaryFieldRef());

            runTime.write();
        }

        runTime.printExecutionTime(Info);
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
