"""

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from fluidfoam import readof as rdf
import os

plt.rcParams["font.size"] = 15


foamTimes = os.popen('foamListTimes').read()
timeList = foamTimes.split('\n')[:-1]
timeArr = np.array([float(t) for t in timeList])

ntimes = len(timeList)

# color list for results
colors = ["cornflowerblue", "tomato"]
cm = LinearSegmentedColormap.from_list(
        "Custom", colors, N=ntimes)
colors = cm(np.arange(0, cm.N))

maxBeta = np.zeros(ntimes)
maxBetaPY = np.zeros(ntimes)


fig, (axZb, axQb, axBeta) = plt.subplots(nrows=3, figsize=(8.3, 9))

for i, time in enumerate(timeList):
    Xbed, Ybed, Zbed = rdf.readmesh(
        "./", time, boundary="bed", verbose=False)
    axZb.plot(Xbed, Zbed, color=colors[i])
    axZb.set_ylabel(r"$z_b\,[m]$")

    # avalanche bedload flux
    qavX, qavY, qavZ = rdf.readvector(
        "./", time, "qavVf", boundary="bed", verbose=False)
    magQav = np.sqrt(qavX**2 + qavZ**2)

    axQb.plot(Xbed, magQav, color=colors[i])
    axQb.set_ylabel(r"$q_{av}\,[m^2.s^{-1}]$")

    betaOF = rdf.readscalar(
        "./", time, "betaVf", boundary="bed", verbose=False)
    maxBeta[i] = np.max(betaOF)

    # compute beta from qav
    gradZb = np.arctan(qavZ/qavX)
    betaPY = np.arctan(np.abs(gradZb))
    maxBetaPY[i] = np.max(betaPY)


axBeta.plot(
    timeArr, maxBeta * 180 / np.pi, color="#D55E00",
    label="OpenFOAM")
axBeta.plot(
    timeArr, maxBetaPY * 180 / np.pi, ls="dashed",
    color="#0072B2", label=r"from $q_{av}$")

axQb.set_xlabel("x [m]")

axBeta.set_ylabel(r"$max(\beta)$ in degrees")
axBeta.set_xlabel(r"time [s]")
axBeta.legend()

fig.tight_layout()

plt.show()
