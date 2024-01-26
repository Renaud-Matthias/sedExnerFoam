"""
Plot velocity and pression colormap for submerged jet case
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

Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase, time_name=time)

Xbed, Ybed, Zbed = rdf.readmesh(
    pathCase,
    time_name=time,
    boundary="sedimentBed"
)

# velocity
U = rdf.readfield(pathCase, time_name=time, name="U")
# pressure
p = rdf.readscalar(pathCase, time_name=time, name="p")


fig, axs = plt.subplots(2, 1, layout="constrained")

Xapron = np.linspace(0, Lapron, 10)
Zapron = np.zeros_like(Xapron)
# y coordinates for bottom of graph
Zbotbed = np.zeros_like(Xbed) - zminplot
Zbotapron =  np.zeros_like(Xapron) - zminplot

# plot p
imp = axs[0].tripcolor(
    Xmesh, Zmesh, p, cmap="Spectral_r", shading="gouraud")
cbarp = fig.colorbar(imp, ax=axs[0])
axs[0].set_ylabel("z [m]")
cbarp.set_label(r"p [$m^2.s^{-2}$]")
axs[0].set_xlim(None, None)
axs[0].set_ylim(None, None)

# plot mag(U)
imU = axs[1].tripcolor(
    Xmesh, Zmesh, np.sqrt(U[0]**2 + U[2]**2), cmap="Spectral_r", shading="gouraud")
cbarU = fig.colorbar(imU, ax=axs[1])
axs[1].set_xlabel("x [m]")
axs[1].set_ylabel("z [m]")
cbarU.set_label(r"U [$m.s^{-1}$]")

for ax in axs:
    ax.fill_between(Xbed, Zbed, Zbotbed, color="peru")
    ax.fill_between(Xapron, Zapron, Zbotapron, color="slategray")
    ax.plot(Xbed, Zbed, color="black")
    ax.plot(Xapron, Zapron, color="black")
    ax.plot([Lapron, Lapron], [-zminplot, 0], color="black")
    ax.set_xlim(min(Xmesh), max(Xmesh))
    ax.set_ylim(-zminplot, max(Zmesh))

plt.show()
