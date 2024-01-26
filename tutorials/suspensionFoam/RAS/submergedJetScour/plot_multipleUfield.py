"""
Plot velocity colormap for submerged jet case for different cases
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from fluidfoam import readof as rdf
import os

# of case and time
pathCase = "./"
time = "latestTime"

# geometrical parameters
# apron length
Lapron = 0.05
# zmin for plots
zminplot = 0.05

timeList = ["3_mpm", "3"]
nt = len(timeList)

fig, axs = plt.subplots(nt, 1, layout="constrained")

for i, t in enumerate(timeList):
    Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase, time_name=t)

    Xbed, Ybed, Zbed = rdf.readmesh(
        pathCase,
        time_name=t,
        boundary="sedimentBed"
    )

    # velocity
    U = rdf.readfield(pathCase, time_name=t, name="U")

    Xapron = np.linspace(0, Lapron, 10)

    Zapron = np.zeros_like(Xapron)
    # y coordinates for bottom of graph

    Zbotbed = np.zeros_like(Xbed) - zminplot
    Zbotapron =  np.zeros_like(Xapron) - zminplot

    # plot mag(U)
    imU = axs[i].tripcolor(
        Xmesh, Zmesh, np.sqrt(U[0]**2 + U[2]**2),
        cmap="Spectral_r", vmax=1.9, shading="gouraud")
    cbarU = fig.colorbar(imU, ax=axs[i])
    axs[i].set_xlabel("x [m]")
    axs[i].set_ylabel("z [m]")
    cbarU.set_label(r"U [$m.s^{-1}$]")
    axs[i].fill_between(Xbed, Zbed, Zbotbed, color="peru")
    axs[i].fill_between(Xapron, Zapron, Zbotapron, color="slategray")
    axs[i].plot(Xbed, Zbed, color="black")
    axs[i].plot(Xapron, Zapron, color="black")
    axs[i].plot([Lapron, Lapron], [-zminplot, 0], color="black")
    axs[i].set_xlim(min(Xmesh), max(Xmesh))
    axs[i].set_ylim(-zminplot, max(Zmesh))

plt.show()
