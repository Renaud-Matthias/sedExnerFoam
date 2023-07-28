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

if (initialBedform=="gaussian")
{
    dimensionedScalar Hdune(dimLength);
    dimensionedScalar x0dune(dimLength);
    dimensionedScalar Sdune(dimLength);

    Hdune.value() = readScalar(bedformDict.lookup("duneHeight"));
    x0dune.value() = readScalar(bedformDict.lookup("dunePos"));
    Sdune.value() = readScalar(bedformDict.lookup("duneWidth"));
    
    forAll(aMesh.areaCentres(), i)
    {
        scalar x = aMesh.areaCentres()[i].component(0);
        scalar X = (x - x0dune.value()) / Sdune.value();
        Zb[i] = Hdune.value() * Foam::exp(-pow(X, 2));
    }
}

else if (initialBedform=="tangente")
{
    dimensionedScalar x0bf(dimLength);
    dimensionedScalar Lbf(dimLength);
    scalar thetabf = 0;

    // position related to initial bed form
    x0bf.value() = readScalar(bedformDict.lookup("bedFormPos"));
    // initial bedorm width
    Lbf.value() = readScalar(bedformDict.lookup("bedFormWidth"));
    // initial maximal slope (in degrees)
    thetabf = readScalar(bedformDict.lookup("angleIni"));

    scalar slopeCI = Foam::tan(Foam::degToRad(thetabf));
    scalar alphaCI = (2/Lbf.value()) * Foam::tan(2*(1-0.9));
    dimensionedScalar heightbf(dimLength);
    heightbf.value() = (2 / alphaCI) * Foam::tan(slopeCI);
    
    forAll(aMesh.areaCentres(), i)
    {
        scalar x = aMesh.areaCentres()[i].component(0);
        Zb[i] = heightbf.value()
            * (1 - 0.5 * Foam::atan(alphaCI* (x - x0bf.value())));
    }
}

else if (initialBedform=="cone")
{
    dimensionedScalar x0cone(dimLength);
    dimensionedScalar Wcone(dimLength);
    scalar thetaCone = 0;

    // x coordinates of cone center
    x0cone.value() = readScalar(bedformDict.lookup("conePos"));
    // hald of cone width
    Wcone.value() = readScalar(bedformDict.lookup("coneWidth"));
    // angle of the slope (in degrees)
    thetaCone = readScalar(bedformDict.lookup("angle"));

    forAll(aMesh.areaCentres(), i)
    {
        scalar x = aMesh.areaCentres()[i].component(0);
        scalar zbinf = (x - x0cone.value() + Wcone.value())
            * Foam::tan(degToRad(thetaCone))
            * pos(x - x0cone.value() + Wcone.value()) * neg(x - x0cone.value());
        scalar zbsup = -(x - x0cone.value() - Wcone.value())
            * Foam::tan(degToRad(thetaCone))
            * pos0(x - x0cone.value())
            * neg(x - x0cone.value() - Wcone.value());
        Zb[i] = zbinf + zbsup;
    }
}
