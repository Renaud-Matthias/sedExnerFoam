"""
1DSedim test
settling of particules with hindrance effects
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from fluidfoam import readof as rdf
import os


plt.rcParams["font.size"] = 15

foamTimes = os.popen("foamListTimes -withZero").read()
timeList = foamTimes.split("\n")[:-1]
ntimes = len(timeList)

# color list for results
colors = ["#0072B2", "#D55E00"]
cm = LinearSegmentedColormap.from_list(
        "Custom", colors, N=ntimes)
colors = cm(np.arange(0, cm.N))


def C0(Y):
    """initial concentration profile"""
    Cout = 0.5 * (1+np.tanh(10*(Y-0.054)/(0.049-Y)))
    return np.where(Y > 0.049, 0.5*Cout, 0.5)


# load results from previous simulation
Ydata, *csData = np.loadtxt(
    "dataSedim.txt", delimiter=";", unpack=True)

Ycells = rdf.readmesh("./")[1]

# get cell volumes
try:
    Vcells = rdf.readscalar(
        "./", "0", "V", verbose=False)
except FileNotFoundError:
    os.system("postProcess -func writeCellVolumes > /dev/null")
    Vcells = rdf.readscalar(
        "./", "0", "V", verbose=False)

# mass of suspended sediments
massSusp = np.zeros(ntimes)

fig, (axCs, axErr) = plt.subplots(ncols=2, figsize=(10, 6))

for i, t in enumerate(timeList):
    if t == "0":
        Cs = C0(Ycells)
    else:
        Cs = rdf.readscalar("./", t, "Cs")
        axCs.plot(
            csData[i-1], Ydata, marker="x", color="black")
        csErr = Cs - csData[i-1]
        axErr.scatter(csErr, Ycells, color=colors[i])
    axCs.scatter(Cs, Ycells, marker="o", color=colors[i])

axCs.set_ylabel("y [m]")
axCs.set_xlabel(r"$c_s$")
axCs.grid()

axErr.set_xlabel("$c_s$ error")
axErr.grid()

fig.tight_layout()

plt.show()
