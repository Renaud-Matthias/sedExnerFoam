"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print("- running test equilibrium channel")

time = "latestTime"

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

shieldsSolv = rdf.readvector("./", time, "shieldsVf",
                             verbose=False, boundary="bed")[0, 0]
qbSolv = rdf.readvector("./", time, "qbVf",
                        verbose=False, boundary="bed")[0, 0]
tauOf = -rdf.readvector("./", time, "wallShearStress",
                        verbose=False, boundary="bed")[0, 0]

shieldsOf = tauOf / ((rhoS/rhoF - 1) * g * dS)

success = True
tol = 1e-4

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
    print("problem, bedload value")

assert success
