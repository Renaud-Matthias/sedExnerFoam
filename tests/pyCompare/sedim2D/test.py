"""
2D sedimentation cases, every time step saved
compute the deposition flux and the associated bed elevation
check if it corresponds to sedExnerFoam results
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- sedim 2D, deposition/bed elevation --- ")

verbose = False

success = True
tolMass = 1e-6

foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
timeArr = np.array([float(t) for t in timeList])

ntimes = len(timeList)

rhoS = 2650.  # sediment density m3/s
CsMax = 0.57  # maximum sediment volume fraction
Hwater = 0.1  # water height
Lx = 0.1  # domain length
Cs0 = 0.05  # initial sediment volume fraction

# initial suspended mass
m0 = rhoS * Cs0 * Lx * Hwater

Xfaces, dummy, Zbed = rdf.readmesh(
    "./", boundary="bed", verbose=False)
nFaces = len(Xfaces)
nEdges = nFaces + 1
Xedges = np.zeros(nEdges)
Xedges[1:-1] = 0.5 * (Xfaces[1:] + Xfaces[:-1])
Xedges[0] = Xfaces[0] - (Xedges[1] - Xfaces[0])
Xedges[-1] = Xfaces[-1] + (Xfaces[-1] - Xedges[-2])
# faces width
dX = Xedges[1:] - Xedges[:-1]

dzEdges = np.zeros(nEdges)

dzbErr = np.zeros(ntimes-1)

ws = rdf.readvector("./", "latestTime", "Ws", verbose=False)[2, 0]
if verbose:
    print(f"settling velocity: {ws*100} cm/s")

# mass increments, from OpenFOAM and computed with python
dmBedOF = np.zeros(ntimes-1)
dmBedPY = np.zeros(ntimes-1)

# initial bed velocity
wBed = np.zeros(nFaces)

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(timeList[i+1]) - float(timeList[i])
    if verbose:
        print(f"\n- time, t = {t} s, dt = {dt} s")
    # read bed pos from mesh a t and t+dt
    # read bed pos from mesh a t and t+dt
    zbFaces_t = rdf.readmesh(
        "./", t, boundary="bed", verbose=False)[2]
    zbFaces_tnext = rdf.readmesh(
        "./", tnext, boundary="bed", verbose=False)[2]
    dzbOF = zbFaces_tnext - zbFaces_t
    # read Cs, sediment volume fraction
    Cs = rdf.readscalar("./", t, "Cs", verbose=False)
    # compute deposition flux
    Dep = (-ws + wBed) * Cs[:nFaces]
    dzbPY = Dep * dt / CsMax
    dzbErr[i] = np.sum(dzbOF - dzbPY)

    # mass in kg/m
    dmBedOF[i] = rhoS * CsMax * np.sum(dzbOF * dX)
    dmBedPY[i] = rhoS * CsMax * np.sum(dzbPY * dX)
    # update bed velocity
    wBed = (zbFaces_tnext - zbFaces_t) / dt


massBedOF = np.cumsum(dmBedOF)
massBedPY = np.cumsum(dmBedPY)

# relative error on mass at each time step
relErr = (massBedOF - massBedPY) / m0
maxRelErr = np.max(np.abs(relErr))

print(f"maximum relative error on mass: {maxRelErr * 100} %")

if maxRelErr > tolMass:
    success = False
    print(
        f"error! relative error on bed mass = "
        + f"{100 * maxRelErr} %")
    print(f"tolerance is {100 * tolMass} %")
else:
    print("mass conservation OK")


assert success
