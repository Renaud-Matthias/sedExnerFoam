"""
compare solution computed with python to sedExnerFoam result
"""

import numpy as np
from fluidfoam import readof as rdf
import matplotlib.pyplot as plt
import os

plt.rcParams["font.size"] = 15


foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
timeArr = np.array([float(t) for t in timeList])

print("time list: ", timeList)

ntimes = len(timeList)

rhoS = 2650.  # sediment density m3/s
CsMax = 0.57  # maximum sediment volume fraction
Hwater = 0.1  # water height
Lx = 0.1  # domain length
Cs0 = 0.05  # initial suspended sediment volume fraction

# initial suspended mass [kg/m]
m0 = rhoS * Cs0 * Lx * Hwater
print(f"initial mass : {m0} kg.m-1")

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
print(f"settling velocity: {ws*100} cm/s")

# mass increments, from OpenFOAM and computed with python
dmBedOF = np.zeros(ntimes-1)
dmBedPY = np.zeros(ntimes-1)

# initial bed velocity
wBed = np.zeros(nFaces)

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(timeList[i+1]) - float(timeList[i])
    print(f"\n- time, t = {t} s, dt = {dt} s")
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


fig, (axM, axErr) = plt.subplots(nrows=2, figsize=(8.3, 7))

fig.suptitle("mass conservation")

axM.plot(timeArr[:-1], massBedOF, label="bed mass OF")
axM.plot(timeArr[:-1], massBedPY, ls="dashed", label="bed mass PY")
axM.axhline(m0, color="black", ls="dotted", label="total mass")

axErr.plot(timeArr[:-1], relErr*100)

axM.legend()
axM.set_ylabel(r"mass $[kg.m{-1}]$")

axErr.set_ylabel("relative mass error in %")
axErr.set_xlabel("time [s]")

fig.tight_layout()
plt.show()
