"""
compare solution computed with python to sedExnerFoam result
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running AdamsBashforth 0 compare slope avalanche --- ")

foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
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

qbEdges = np.zeros(nEdges)
dzEdges = np.zeros(nEdges)

success = True
tol = 0.2  # tolerance in %

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
    qb_x = rdf.readvector(
        "./", t, "qbVf", boundary="bed", verbose=False)[0]
    # qb on edges, linear interpolation
    qbEdges[1:-1] = 0.5 * (qb_x[1:] + qb_x[:-1])
    # zero flux boundary condition for qb
    qbEdges[0], qbEdges[-1] = 0, 0
    # bed level increment between t and t+dt, exner equation
    dzb = (1/CsMax) * (qbEdges[:-1] - qbEdges[1:]) * dt / dX
    # bed level increment on edges, linear interpolation
    dzEdges[1:-1] = 0.5 * (dzb[1:] + dzb[:-1])
    dzEdges[0], dzEdges[-1] = dzb[0], dzb[-1]
    # bed level increment on faces, linear interpolation
    dzb = 0.5 * (dzEdges[1:] + dzEdges[:-1])
    Zbed += dzb

    zBedErr = 100 * np.abs(dzbOF - dzb)/np.max(np.abs(dzb))
    maxErr = np.max(zBedErr)
    print(f"max relative error on bed increment = {round(maxErr, 5)} %")

    if (maxErr > tol):
        success = False

if success:
    print("test passed")
assert success
