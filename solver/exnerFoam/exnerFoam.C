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
        \ddt{Zb} + \div(phis)
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
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Exner equation solver."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"
    #include "createFaFields.H"

    Info << "number of areas" << aMesh.nFaces() << endl;

    // set initial condition for dune tutorial
    // should later be directly in the tutorial case
    forAll(aMesh.areaCentres(), i)
    {
        scalar x = aMesh.areaCentres()[i].component(0);
        Zb[i] = 0.2*Foam::exp(-0.1*pow(x-10, 2));
    }
    
    #include "createVolFields.H"

    dimensionedVector Q("Q", dimVelocity*dimLength, vector(1, 0, 0));
    dimensionedScalar H("H", dimLength, 1);

    Info << "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.value() << endl;

        forAll(aMesh.areaCentres(), i)
        {
            scalar qb = Q.value().x()/(H.value()-Zb[i]); // explicit
            //scalar qb = Q.value().x()/Zb[i]*(H.value()-Zb[i]); //implicit
            Qb[i] = vector(qb, 0, 0);
        }

        faScalarMatrix ZbEqn
            (
                fam::ddt(Zb)
                + fac::div(Qb) //explicit
                //+ fam::div(phib,Zb) //implicit
            );

        ZbEqn.solve();

        //Zb.correctBoundaryConditions();
        //Qb.correctBoundaryConditions();

        Info<< "bed elevation = "
            << Zb.weightedAverage(aMesh.S()).value()
            << "  Min(Zb) = " << gMin(Zb)
            << "  Max(Zb) = " << gMax(Zb)
            << endl;

        if (runTime.writeTime())
        {
            vsm.mapToVolume(Zb, Zbvf.boundaryFieldRef());

            runTime.write();
        }

        runTime.printExecutionTime(Info);
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
