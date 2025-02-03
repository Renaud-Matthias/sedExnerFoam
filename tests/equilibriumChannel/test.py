"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print(" --- running test equilibrium channel --- ")

time = "50"

dS = 0.2e-3  # sediment diameter
rhoS = 2650.  # sediment density
rhoF = 1000.  # fluid density
g = 9.81  # gravity acceleration
critShields = 0.047


def bedloadMPM(shields):
    """
    Return bedload computed using Meyer-Peter Muller formula
    """
    einNum = np.sqrt((rhoS/rhoF - 1) * g * dS**3)
    return 8 * einNum * (shields - critShields)**1.5


Zmesh = rdf.readmesh("./", verbose=False)[2]

Cs = rdf.readscalar("./", time, "Cs", verbose=False)

shieldsSolv = rdf.readvector(
    "./", time, "shieldsVf", verbose=False, boundary="bed")[0, 0]
qbSolv = rdf.readvector(
    "./", time, "qbVf", verbose=False, boundary="bed")[0, 0]
tauOf = -rdf.readvector(
    "./", time, "wallShearStress", verbose=False, boundary="bed")[0, 0]

shieldsOf = tauOf / ((rhoS/rhoF - 1) * g * dS)

success = True
tol = 1e-5

err = np.abs((shieldsSolv - shieldsOf) / shieldsSolv)
# test shields value
if err > tol:
    success = False
    print("problem in shields number value")
# test qb value with Meyer-Peter law
qbTest = bedloadMPM(shieldsOf)
err = np.abs((qbSolv - qbTest) / qbSolv)
if err > tol and success is True:
    success = False
    print("error, bedload value != Meyer-Peter formula")


zDatat, CsData = np.loadtxt("./dataCs_50s.txt", unpack=True, delimiter=";")
for csim, cdat in zip(Cs, CsData):
    err = np.abs(csim - cdat)
    if err > tol:
        success = False
        print("error, new results differs from old simulation")

assert success
