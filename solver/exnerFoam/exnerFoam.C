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
#include "faCFD.H"
#include "unitConversion.H"

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
    #include "createFaFields.H"

    Info << "number of areas" << aMesh.nFaces() << endl;

    // set initial condition for dune tutorial
    // initiate Zb level, dune or triangle for slide
    if (setupType=="dune")
    {
        forAll(aMesh.areaCentres(), i)
        {
            scalar x = aMesh.areaCentres()[i].component(0);
            scalar X = (x - x0dune.value()) / Sdune.value();
            //scalar y = aMesh.areaCentres()[i].component(1);
            //scalar distCenterSqr = pow(x - 2, 2) + pow(y, 2);
            //Zb[i] = 0.1 * Foam::exp(-distCenterSqr / 0.36);
            Zb[i] = Hdune.value() * Foam::exp(-pow(X, 2));
        }
    }

    else if (setupType=="avalanche")
    {
        scalar slopeCI = Foam::tan(Foam::degToRad(thetabf));
        
        scalar alphaCI = (2/Lbf.value()) * Foam::tan(2*(1-0.9));
        
        dimensionedScalar heightbf(dimLength);
        heightbf.value() = (2 / alphaCI) * Foam::tan(slopeCI);

        forAll(aMesh.areaCentres(), i)
        {
            scalar x = aMesh.areaCentres()[i].component(0);
            /*
            scalar zl = slopeCI * (x - x0cone.value() + coneW.value())
                * pos(x - x0cone.value() + coneW.value())
                * neg0(x - x0cone.value());
            scalar zr = -slopeCI * (x - x0cone.value() - coneW.value())
                * pos(x - x0cone.value())
                * neg(x - x0cone.value() - coneW.value());
            Zb[i] = zl + zr;*/
            Zb[i] = heightbf.value()
                * (1 - 0.5 * Foam::atan(alphaCI* (x - x0bf.value()))); 
        }
    }
    
    #include "createVolFields.H"

    Info << "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.value() << endl;

        theta = Foam::tan(Foam::mag(fac::grad(Zb)));

        if (setupType=="dune")
        {
            forAll(aMesh.areaCentres(), i)
            {
                // explicit
                scalar qb = alpha *
                    Foam::pow(
                        Qwater.value().x() / (Hwater.value() - Zb[i]), beta);
                // implicit
                //scalar qb = Q.value().x()/Zb[i]*(H.value()-Zb[i]);
                Qb[i] = vector(qb, 0, 0);
            }
        }
        else if (setupType=="avalanche")
        {
            forAll(aMesh.areaCentres(), i)
            {
                /*scalar thetaLocal = theta[i] - Foam::degToRad(thetaRep);
                scalar da = Foam::tan(
                    Foam::mag(thetaLocal)) * pos(thetaLocal);
                Da[i] = da;*/
                if (theta[i] > Foam::degToRad(thetaRep))
                {
                    Da[i] = 0.1;
                }
                else
                {
                    Da[i] = 0.0;
                }
            }
        }

        faScalarMatrix ZbEqn
            (
                fam::ddt(Zb)
                + fac::div(Qb) //explicit
                - fac::laplacian(Da, Zb)
                //+ fam::div(phib,Zb) //implicit
            );

        ZbEqn.solve();

        Zb.correctBoundaryConditions();
        Qb.correctBoundaryConditions();

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
