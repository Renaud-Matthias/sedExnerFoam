"""
compare solution computed with python to sedExnerFoam result
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running AdamsBashforth 0 compare slope avalanche --- ")

verbose = False

foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
if verbose:
    print("time list: ", timeList)

ntimes = len(timeList)

rhoS = 2650.  # sediment density m3/s
CsMax = 0.57  # maximum sediment volume fraction

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

qavEdges = np.zeros(nEdges)
dzEdges = np.zeros(nEdges)

success = True
tol = 0.2  # tolerance in %

for i in range(ntimes-1):
    t, tnext = timeList[i], timeList[i+1]
    dt = float(timeList[i+1]) - float(timeList[i])
    if verbose:
        print(f"\n- time, t = {t} s, dt = {dt} s")
    # read bed pos from mesh a t and t+dt
    zbFaces_t = rdf.readmesh(
        "./", t, boundary="bed", verbose=False)[2]
    zbFaces_tnext = rdf.readmesh(
        "./", tnext, boundary="bed", verbose=False)[2]
    dzbOF = zbFaces_tnext - zbFaces_t
    qav_x = rdf.readvector(
        "./", t, "qavVf", boundary="bed", verbose=False)[0]
    # qb on edges, linear interpolation
    qavEdges[1:-1] = 0.5 * (qav_x[1:] + qav_x[:-1])
    # zero flux boundary condition for qb
    qavEdges[0], qavEdges[-1] = 0, 0
    # bed level increment between t and t+dt, exner equation
    dzb = (1/CsMax) * (qavEdges[:-1] - qavEdges[1:]) * dt / dX
    # bed level increment on edges, linear interpolation
    dzEdges[1:-1] = 0.5 * (dzb[1:] + dzb[:-1])
    dzEdges[0], dzEdges[-1] = dzb[0], dzb[-1]
    # bed level increment on faces, linear interpolation
    dzb = 0.5 * (dzEdges[1:] + dzEdges[:-1])
    Zbed += dzb

    zBedErr = 100 * np.abs(dzbOF - dzb)/np.max(np.abs(dzb))
    maxErr = np.max(zBedErr)
    if verbose:
        print(
            f"max relative error on bed increment = {round(maxErr, 5)} %")
    if (maxErr > tol):
        success = False
        print(
            f"error! relative error on bed level increment = "
            + f"{100*maxxErr} %")
        print(f"tolerance is {100 * tol} %")

assert success
