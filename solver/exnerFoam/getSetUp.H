/*---------------------------------------------------------------------------* \
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

Info << "Reading setupProperties" << endl;

IOdictionary setupProperties
(
    IOobject
    (
        "setupProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// initial bedform
dictionary bedformDict(setupProperties.subDict("bedform"));
word initialBedform(bedformDict.lookup("type"));

word bedload(setupProperties.lookup("bedload"));
word avalanche(setupProperties.lookup("avalanche"));
//word setupType(setupProperties.lookup("setupType"));

// variables for dune transport case
scalar alpha = 0;
scalar beta = 0;
dimensionedScalar Hwater(dimLength);
dimensionedVector Qwater(dimensionSet(0, 2, -1, 0, 0, 0, 0));

// variables for avalanche case
word avalancheModel;
dimensionedScalar DaMax(dimensionSet(0, 2, -1, 0, 0, 0, 0));
scalar thetaRep = 0;  // angle of repose (in degrees)
dimensionedScalar x0bf(dimLength);
dimensionedScalar Lbf(dimLength);

if (bedload=="on")
{
    Info << "Select dune transport case" << endl;
    Hwater.value() = readScalar(setupProperties.lookup("Hwater"));
    Qwater.value() = vector(setupProperties.lookup("Qwater"));
    alpha = readScalar(setupProperties.lookup("alpha"));
    beta = readScalar(setupProperties.lookup("beta"));
}

if (avalanche=="on")
{
    Info << "Select avalanche case" << endl;
    avalancheModel = word(setupProperties.lookup("avalancheModel"));
    DaMax.value() = readScalar(setupProperties.lookup("DaMax"));
    thetaRep = readScalar(setupProperties.lookup("reposeAngle"));
}

//- resolution of exner equation
// choice between normal openfoam scheme (implicit/explicit solve)
// or Adam Bashforth scheme (explicit)
word eqRes;
if (setupProperties.found("resolution"))
{
    eqRes = word(setupProperties.lookup("resolution"));
}
else
{
    eqRes = "explicit";
}

if (eqRes!="explicit" and eqRes!="implicit" and eqRes!="custom")
{
    FatalError << "wrong keyword in setupProperties dictionary! "
        << "resolution entry, possible values are : "
        << "explicit implicit custom" << endl;
    Info << abort(FatalError) << endl;
}

word schemeCustom;
if (eqRes=="custom")
{
    schemeCustom = word(setupProperties.lookup("schemeCustom"));
}
