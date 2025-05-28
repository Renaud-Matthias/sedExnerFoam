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

rhoS = 2600.  # sediment density m3/s
CsMax = 0.6  # maximum sediment volume fraction
Hwater = 0.5  # water column height
Cs0 = 0.2  # initial sediment volume fraction

# initial suspended mass (kg.m-2)
m0 = rhoS * Cs0 * Hwater

ws = rdf.readvector("./", "latestTime", "Ws", verbose=False)[2, 0]
print(f"settling velocity: {ws*100} cm/s")

massBedOF = np.zeros(ntimes-1)
massBedPY = np.zeros(ntimes-1)

# initial bed velocity
wBed = 0.

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(timeList[i+1]) - float(timeList[i])
    print(f"\n- time, t = {t} s, dt = {dt} s")
    # read bed pos from mesh a t and t+dt
    zbFaces_t = rdf.readmesh(
        "./", t, boundary="bed", verbose=False)[2][0]
    zbFaces_tnext = rdf.readmesh(
        "./", tnext, boundary="bed", verbose=False)[2][0]
    dzbOF = zbFaces_tnext - zbFaces_t
    # read Cs, sediment volume fraction
    Cs = rdf.readscalar("./", t, "Cs")[0]
    # compute deposition
    Dep = (-ws + wBed) * Cs
    # mass increment deposition on bed
    dm = rhoS * Dep * dt

    massBedOF[i] = rhoS * CsMax * zbFaces_tnext
    massBedPY[i] = dm
    # update bed velocity
    wBed = (zbFaces_tnext - zbFaces_t) / dt

massBedPY = np.cumsum(massBedPY)

relErr = (massBedOF - massBedPY) / m0
maxRelErr = np.max(np.abs(relErr))

fig, (axM, axErr) = plt.subplots(nrows=2)

axM.plot(timeArr[:-1], massBedOF, color="#0072B2", label="OF")
axM.plot(timeArr[:-1], massBedPY, ls="dashed", color="#D55E00", label="PY")

axErr.plot(timeArr[:-1], 100*relErr)

axM.legend()
axM.set_ylabel(r"mass $[kg.m^{-2}]$")

axErr.set_ylabel("relative error in %")
axErr.set_xlabel("time [s]")

fig.tight_layout()
plt.show()
