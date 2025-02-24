"""
1D simulation of erosion and deposition
test if sediment mass is conserved
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running 1D erosion test --- ")

success = True
tolMass = 1e-8

foamTimes = os.popen("foamListTimes -withZero").read()
timeList = foamTimes.split("\n")[:-1]
timeArr = np.array([float(t) for t in timeList])
ntimes = len(timeList)

rhoS = 2650.  # sediment density
CsMax = 0.57  # maximum sediment volume fraction

# domain and mesh dimensions
domHeight = 1.
dx, dy = 0.01, 0.001

bedArea = dx * dy

Msusp = np.zeros(ntimes)  # suspended mass
Mbed = np.zeros(ntimes)  # bed mass

Zmesh = rdf.readmesh("./", "latestTime", verbose=False)[2]
nCells = len(Zmesh)

for i, t in enumerate(timeList):
    zb = rdf.readmesh("./", t, boundary="bed", verbose=False)[2][0]
    Mbed[i] = CsMax * rhoS * zb * dx * dy
    CsField = rdf.readscalar("./", t, "Cs", verbose=False)
    if CsField.shape == (1,):
        CsField = np.ones(nCells) * CsField
    try:
        Vcells = rdf.readscalar(
            "./", t, "V", verbose=False)
    except FileNotFoundError:
        os.system("postProcess -func writeCellVolumes > /dev/null")
        Vcells = rdf.readscalar(
            "./", t, "V", verbose=False)
    Msusp[i] = rhoS * np.sum(CsField * Vcells)

Mtotal = Msusp + Mbed
errMass = np.abs(Mtotal) / np.max(Msusp)
endSimuErr = errMass[-1]

# test conservation of mass
if endSimuErr > tolMass:
    success = False
    print(f"error! maximum relative error on mass: {100*endSimuErr} %")
    print(f"tolerance is {100*tolMass} %")
else:
    print("mass conservation OK")

assert success
