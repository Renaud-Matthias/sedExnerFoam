"""
Check mass conservation, between suspended load and bed during deposition
"""

import numpy as np
import os
import matplotlib.pyplot as plt

import fluidfoam.readof as rdf


def myReadmesh(path, t):
    Zcells = rdf.readmesh(path, time_name=t, verbose=False)[2]
    ncells = len(Zcells)
    zbed = rdf.readmesh(
        path, time_name=t, boundary="bed", verbose=False)[2][0]
    ztop = rdf.readmesh(
        path, time_name=t, boundary="top", verbose=False)[2][0]
    Zfaces = np.zeros(ncells+1)
    print(zbed, ztop)
    Zfaces[0], Zfaces[-1] = zbed, ztop
    for i in range(1, ncells):
        Zfaces[i] = 1.5 * Zcells[i-1] - 0.5 * Zfaces[i-1]
    return Zcells, Zfaces


foamTimes = os.popen('foamListTimes -withZero').read()

timeList = foamTimes.split('\n')[:-1]

ntimes = len(timeList)

rhoS = 2600  # sediment density
porosity = 0.6

# domain and mesh dimensions
domHeight = 1.
domWidth = 0.1
ncells = len(rdf.readmesh("./", time_name="0")[0])

bedArea = domWidth**2

print("-domain dimensions")
print(f"height {domHeight} m")
print(f"area {round(bedArea, 5)} m2")
print(f"number of cells {ncells}")

massSuspension = np.zeros(ntimes)
massBed = np.zeros(ntimes)
totMass = np.zeros(ntimes)
relMassErr = np.zeros(ntimes)

for i, t in enumerate(timeList):
    Zcells, Zfaces = myReadmesh("./", t)
    print(f"bed elevation: {Zfaces[0]} m")
    if t == "0":
        Cs = np.ones(ncells) * rdf.readscalar(
            "./", time_name=t, name="Cs", verbose=False)[0]
    else:
        Cs = rdf.readscalar("./", time_name=t, name="Cs")
    try:
        cellVolumes = rdf.readscalar("./", time_name=t, name="V")
    except:
        os.system("postProcess -func writeCellVolumes")
        cellVolumes = rdf.readscalar("./", time_name=t, name="V")
    # bed is only one face, mesh is 1D
    massBed[i] = rhoS * porosity * Zfaces[0] * bedArea
    massSuspension[i] = rhoS * np.sum(Cs * cellVolumes)
    totMass[i] = massBed[i] + massSuspension[i]
    relMassErr[i] = 100 * (totMass[i] - totMass[0]) / totMass[0]

# total mass of sediments
print(f"initial mass: {round(totMass[0], 5)} kg, "
      + f"final mass: {round(totMass[-1], 5)} kg")

fig = plt.figure(figsize=(7, 6), layout="constrained")
ax0 = fig.add_subplot(2, 1, 1)
ax0.plot(
    timeList, totMass, ls="solid", lw=1.8, color="#0072B2", label="total")
ax0.plot(
    timeList, massSuspension, ls="dashed", lw=2,
    color="#009E71", label="suspension")
ax0.plot(
    timeList, massBed, ls="dashdot", lw=2, color="#D55E00", label="bed")
ax0.set_ylabel("mass [kg]")
ax0.set_xlabel("t [s]")
ax0.legend()
ax0.grid()

ax1 = fig.add_subplot(2, 1, 2)
ax1.plot(timeList, relMassErr, color="#0072B2")
ax1.set_ylabel("relative mass change, in %")
ax1.set_xlabel("time [s]")
ax1.grid()
plt.show()
