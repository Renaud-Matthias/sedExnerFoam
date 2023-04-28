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

//string setupType("setupType", setupProperties);
word setupType(setupProperties.lookup("setupType"));

// variables for dune transport case
scalar alpha = 0;
scalar beta = 0;
dimensionedScalar Hwater(dimLength);
dimensionedVector Qwater(dimensionSet(0, 2, -1, 0, 0, 0, 0));
dimensionedScalar Hdune(dimLength);
dimensionedScalar x0dune(dimLength);
dimensionedScalar Sdune(dimLength);

// variables for avalanche case
scalar thetaRep = 0;  // angle of repose (in degrees)
scalar thetabf = 0;  // initial angle of cone (in degrees)
dimensionedScalar x0bf(dimLength);  // position related to initial bed form
dimensionedScalar Lbf(dimLength);  // initial bedorm width

if (setupType=="dune")
{
    Info << "Select dune transport case" << endl;
    Hwater.value() = readScalar(setupProperties.lookup("Hwater"));
    Qwater.value() = vector(setupProperties.lookup("Qwater"));
    Info << "flag 1" << endl;
    Hdune.value() = readScalar(setupProperties.lookup("duneHeight"));
    Info << "flag 2" << endl;
    x0dune.value() = readScalar(setupProperties.lookup("dunePos"));
    Sdune.value() = readScalar(setupProperties.lookup("duneWidth"));
    alpha = readScalar(setupProperties.lookup("alpha"));
    beta = readScalar(setupProperties.lookup("beta"));
}

else if (setupType=="avalanche")
{
    Info << "Select avalanche case" << endl;
    word zb0type(setupProperties.lookup("initialBedForm"));
    x0bf.value() = readScalar(setupProperties.lookup("bedFormPos"));
    Lbf.value() = readScalar(setupProperties.lookup("bedFormWidth"));
    thetaRep = readScalar(setupProperties.lookup("reposeAngle"));
    thetabf = readScalar(setupProperties.lookup("angleIni"));
}

else
{
    Info << "in setupProperties :" << endl;
    Info << setupType << " not recognised" << endl;
    //error wordNotFound("false setupType specified");
    //wordNotFound << exit();
}
