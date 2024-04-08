"""
Test for bedload value
"""

import numpy as np
from fluidfoam import readof as rdf

print("- running sedimentation mass conservation")

rhoS = 2650.  # sediment density
porosity = 0.4


def myReadMesh(path, t):
    Zcells = rdf.readmesh(path, time_name=t, verbose=False)[2]
    ncells = len(Zcells)
    zbed = rdf.readmesh(
        path, time_name=t, boundary="bed", verbose=False)[2][0]
    ztop = rdf.readmesh(
        path, time_name=t, boundary="top", verbose=False)[2][0]
    Zfaces = np.zeros(ncells+1)
    Zfaces[0], Zfaces[-1] = zbed, ztop
    for i in range(1, ncells):
        Zfaces[i] = 1.5 * Zcells[i-1] - 0.5 * Zfaces[i-1]
    return Zcells, Zfaces


Zc0, Zf0 = myReadMesh("./", "0")
Zct, Zft = myReadMesh("./", "latestTime")
Cs0 = rdf.readscalar("./", "0", "Cs", verbose=False)
Cst = rdf.readscalar("./", "latestTime", "Cs", verbose=False)

dS = 0.01  # bed face area

initMass = rhoS * dS * np.sum(Cs0 * (Zf0[1:] - Zf0[:-1]))
endMass = rhoS * dS * (
    np.sum(Cst * (Zft[1:] - Zft[:-1]))
    + (1-porosity) * Zft[0])

success = True
tol = 0.01

relMassErr = (endMass - initMass) / initMass
# test shields value
if np.abs(relMassErr) > tol:
    success = False
    print(f"! relative mass change : {relMassErr*100} %")
    print(f"tolerance is {tol*100} %")

assert success
