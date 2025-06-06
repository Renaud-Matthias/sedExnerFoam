"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print(" --- running slope avalanche --- ")

success = True
tolBetaRep = 1e-3

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

assert success
