"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print(" --- running test equilibrium channel --- ")

success = True

tolSh = 1e-5
tolQb = 1e-5
tolCs = 1e-5

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

# Shields number from sedExnerFoam
shieldsSolv = rdf.readvector(
    "./", time, "shieldsVf", verbose=False, boundary="bed")[0, 0]
qbSolv = rdf.readvector(
    "./", time, "qbVf", verbose=False, boundary="bed")[0, 0]
tauOf = -rdf.readvector(
    "./", time, "wallShearStress", verbose=False, boundary="bed")[0, 0]
# Shields number from wallShearStress function
shieldsOf = tauOf / ((rhoS/rhoF - 1) * g * dS)

errSh = np.abs((shieldsSolv - shieldsOf) / shieldsSolv)
# test shields value
if errSh > tolSh:
    success = False
    print(
        "ERROR! Shields number values from sedExnerFoam and "
        + "wallShearStress utility are not matching"
        + f"\nrelative error on Sields number: {errSh}")
else:
    print(f"Shields number value OK, relative error {errSh}")

# test qb value with Meyer-Peter law
qbTest = bedloadMPM(shieldsOf)
errQb = np.abs((qbSolv - qbTest) / qbSolv)
if errQb > tolQb:
    success = False
    print("error, sedExnerFoam bedload value != Meyer-Peter formula")
else:
    print(f"bedload value OK, relative error {errQb}")

# load results from prevous simulation
zDatat, CsData = np.loadtxt("./dataCs_50s.txt", unpack=True, delimiter=";")
errCs = np.abs(Cs - CsData) / CsData
if np.any(errCs > tolCs):
    success = False
    print("ERROR! Cs profile differs from previous reults")
else:
    print(f"Cs profile OK, maximum relative error {np.max(errCs)}")

assert success
