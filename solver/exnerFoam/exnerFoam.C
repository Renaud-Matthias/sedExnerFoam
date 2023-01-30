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
    surfactantFoam

Group
    grpFiniteAreaSolvers

Description
    finiteArea equation solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{Zb} + \div(phis,Zb)
    \f]

    Where:
    \vartable
        Cs      | Passive scalar
        Ds      | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        Cs      | Passive scalar
        Us      | Velocity [m/s]
    \endplaintable

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

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
    #include "createVolFields.H"

    Info << "number of points" << aMesh.nPoints() << endl;

    Info << "area centers " << aMesh.areaCentres().component(0) << endl;
    
    dimensionedScalar Zbmax(dimLength, 0.2);

    dimensionedScalar xMid(dimLength, 5);

    dimensionedScalar duneWidth(dimLength, 0.5);
    
    Zbf = Zbmax*exp(pow((aMesh.areaCentres().component(0)-xMid)/duneWidth, 2));

    Info << "Zb field initiated" << endl;

    Zbf.correctBoundaryConditions();

    Info << "Zb boundary conditions corrected" << endl;
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.value() << endl;

        faScalarMatrix ZbEqn
        (
            fam::ddt(Zbf)
          + fam::div(phis, Zbf)
        );

        ZbEqn.solve();

        if (runTime.writeTime())
        {
            vsm.mapToVolume(Zbf, Zb.boundaryFieldRef());

            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
