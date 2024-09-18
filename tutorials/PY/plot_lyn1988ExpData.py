"""
Plot data from Lyn 1988
"""

import numpy as np
import matplotlib.pyplot as plt
from readLynExpData import readLynParams

pathData = "../DATA/dataLyn1988/"

# water depth for each case
caseList = readLynParams(pathData + "parametersExpLyn.txt")
for case in caseList:
    if case["name"] == "1565":
        case["marker"] = "+"
        case["color"] = "black"
    if case["name"] == "1965":
        case["marker"] = "^"
        case["color"] = "black"
    if case["name"] == "2565":
        case["marker"] = "x"
        case["color"] = "black"
    if case["name"] == "1957":
        case["marker"] = "s"
        case["color"] = "black"


fig, axs = plt.subplots(3, 1, figsize=(8, 8), layout="constrained")

for case in caseList:
    pathC = pathData + "Ceq" + case["name"] + ".txt"
    pathUuf = pathData + "Uuf" + case["name"] + ".txt"
    pathRstress = pathData + "Rstress" + case["name"] + ".txt"
    Zc, C = np.loadtxt(pathC, delimiter=";", unpack=True)
    Zu, Uuf = np.loadtxt(pathUuf, delimiter=";", unpack=True)
    axs[0].scatter(Uuf, Zu, marker=case["marker"],
                   color=case["color"], label="exp " + case["name"])
    axs[1].scatter(C, Zc, marker=case["marker"], color=case["color"])
    if case["name"] != "1565":
        ZRstress, Rstress = np.loadtxt(
            pathRstress, delimiter=";", unpack=True)
        axs[2].scatter(Rstress, ZRstress,
                       marker=case["marker"], color=case["color"])

axs[0].set_ylabel(r"$z^+$")
axs[0].set_xlabel(r"$u/u_f$")
axs[0].set_yscale("log")
axs[0].legend()
axs[0].grid()

axs[1].set_ylabel("z/H")
axs[1].set_xlabel("c")
axs[1].set_xscale("log")
axs[1].set_yscale("log")
axs[1].grid()
axs[1].set_ylim(0.01, 1)
axs[1].set_xlim(1e-7, 0.01)

axs[2].grid()
axs[2].set_ylabel("z/H")
axs[2].set_xlabel(r"$\frac{-u'w'}{u_f^2}$")
axs[2].set_xlim(0, 1)
axs[2].set_ylim(0, 1)

plt.show()
