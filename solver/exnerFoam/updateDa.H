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

if (avalancheModel == "threshold")
{
    forAll(aMesh.areaCentres(), i)
    {
        if (theta[i] > Foam::degToRad(thetaRep))
        {
            Da[i] = DaMax.value();
        }
        else
        {
            Da[i] = 0.0;
        }
    }
}
else if (avalancheModel == "Vinent")
{
    forAll(aMesh.areaCentres(), i)
    {
        Da[i] = DaMax.value() * Foam::pos(theta[i] - Foam::degToRad(thetaRep))
            * (Foam::tanh(Foam::tan(theta[i])
            - Foam::tanh(Foam::tan(Foam::degToRad(thetaRep)))))
            / Foam::tan(theta[i]);
    }
}
else
{
    forAll(aMesh.areaCentres(), i)
    {
        Da[i] = 0.0;
    }
}
