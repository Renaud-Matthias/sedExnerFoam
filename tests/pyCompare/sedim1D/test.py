"""
1D sedimentation cases, every time step saved
compute the deposition flux and the associated bed elevation
check if it corresponds to sedExnerFoam results
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- sedim 1D, deposition/bed elevation --- ")

verbose = False

success = True
tolMass = 5e-6

foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
if verbose:
    print("time list: ", timeList)

ntimes = len(timeList)

rhoS = 2600.  # sediment density m3/s
CsMax = 0.6  # maximum sediment volume fraction
Hwater = 0.5  # water column height
Cs0 = 0.2  # initial sediment volume fraction

# initial suspended mass (kg.m-2)
m0 = rhoS * Cs0 * Hwater

ws = rdf.readvector("./", "latestTime", "Ws", verbose=False)[2, 0]
if verbose:
    print(f"settling velocity: {ws*100} cm/s")

massBedOF = np.zeros(ntimes-1)
massBedPY = np.zeros(ntimes-1)

# initial bed velocity
wBed = 0.

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(timeList[i+1]) - float(timeList[i])
    if verbose:
        print(f"\n- time, t = {t} s, dt = {dt} s")
    # read bed pos from mesh a t and t+dt
    zbFaces_t = rdf.readmesh(
        "./", t, boundary="bed", verbose=False)[2][0]
    zbFaces_tnext = rdf.readmesh(
        "./", tnext, boundary="bed", verbose=False)[2][0]
    dzbOF = zbFaces_tnext - zbFaces_t
    # read Cs, sediment volume fraction
    Cs = rdf.readscalar("./", t, "Cs", verbose=False)[0]
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
if verbose:
    print(f"maximum relative error on mass: {maxRelErr*100} %")

if maxRelErr > tolMass:
    success = False
    print(
        f"error! relative error on bed mass = "
        + f"{100 * maxRelErr} %")
    print(f"tolerance is {100 * tolMass} %")
else:
    print("mass conservation OK")

assert success
