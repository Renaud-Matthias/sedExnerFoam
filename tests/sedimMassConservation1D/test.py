"""
sedimentation of sediment in 1D
test if sediment mass is conserved
all time steps must be saved
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running sedimentation mass conservation --- ")

rhoS = 2600.  # sediment density
CsMax = 0.6  # maximum sediment volume fraction

# domain and mesh dimensions
domHeight = 1.
domWidth = 0.1

bedArea = domWidth**2

ncells = len(rdf.readmesh("./", verbose=False)[2])

foamTimes = os.popen('foamListTimes -withZero').read()
timeList = foamTimes.split('\n')[:-1]
ntimes = len(timeList)


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


massSuspension = np.zeros(ntimes-1)
massBed = np.zeros(ntimes-1)

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(tnext) - float(t)
    if i==0:
        zbed_t = 0.  # initial bed position
    else:
        zbed_t= Zfaces[0]
    
    # get suspended mass
    Cs = rdf.readscalar("./", time_name=t, name="Cs", verbose=False)
    if Cs.shape==(1,):
        Cs = np.ones(ncells) * Cs
    try:
        cellVolumes = rdf.readscalar("./", t, "V")
    except:
        os.system("postProcess -func writeCellVolumes")
        cellVolumes = rdf.readscalar("./", t, "V")
    if cellVolumes.shape==(1,):
        cellVolumes = np.ones(ncells) * cellVolumes
    massSuspension[i] = rhoS * np.sum(Cs * cellVolumes)
    # get bed mass
    Zcells, Zfaces = myReadMesh("./", tnext)
        
    # bed is only one face, mesh is 1D
    massBed[i] = rhoS * CsMax * Zfaces[0] * bedArea


totMass = massBed + massSuspension
# relative mass gain or loss
relMassErr = (totMass - totMass[0]) / totMass[0]

maxRelErr = np.max(np.abs(relMassErr))

success = True
tol = 1e-5
# test shields value
if np.any(np.abs(relMassErr) > tol):
    success = False
    print(f"warning! maximum relative error on mass : {maxRelErr}")
    print(f"tolerance is {tol}")

assert success
