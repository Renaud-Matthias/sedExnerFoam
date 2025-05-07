"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print(" --- running slope avalanche --- ")

success = True
tolBetaRep = 1e-3
tolQav = 1e-4

time = "latestTime"

dS = 0.5e-3  # sediment diameter
rhoS = 2650.  # sediment density
rhoF = 1000.  # fluid density
g = 9.81  # gravity acceleration
repAngle = 32.  # repose angle of granular material (in rad)

# bed position
Xbed, Ybed, Zbed = rdf.readmesh(
    "./", "latestTime", boundary="bed", verbose=False)

tanBed = np.abs((Zbed[2:] - Zbed[:-2]) / (Xbed[2:] - Xbed[:-2]))

maxAngle = (180/np.pi) * np.max(np.arctan(tanBed))

relAngleError = (maxAngle - repAngle) / repAngle
if relAngleError > tolBetaRep:
    success = False
    print(
        f"error! relative difference between repose angle "
        + f"({repAngle} deg) and maximum slope angle ({maxAngle} deg) "
        + f"is {relAngleError*100} %\ntolerance is {tolBetaRep*100} %")
else:
    print("bed slope OK")

# test value of avalanche
Qav0 = 5e-3
betaRep = repAngle * np.pi / 180
qavX, qavY, qavZ = rdf.readvector(
    "./", "0.1", "qavVf", boundary="bed", verbose=False)
magQav = np.sqrt(qavX**2 + qavZ**2)
beta = np.arctan(np.abs(qavZ/qavX))
qavVinent = Qav0 * (np.tanh(np.tan(beta)) - np.tanh(np.tan(betaRep)))
qavVinent /= (1 - np.tanh(np.tan(betaRep)))

relErrQav = np.max(np.abs(magQav - qavVinent) / np.max(qavVinent))
if relErrQav > tolQav:
    success = False
    print(
        f"error! relative error on avalanche bedload flux value: "
        + f"{relErrQav*100} %\ntolerance is {tolQav*100} %")
else:
    print("avalanche flux OK")


assert success
