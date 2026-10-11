import numpy as np
from fluidfoam import readof as rdf
import os

rhoS = 2650.0
CsMax = 0.57
dy = 0.01

foamTimes = os.popen("foamListTimes -withZero").read()
timeList = foamTimes.splitlines()

# Géométrie du fond
Xfaces, _, _ = rdf.readmesh("./", timeList[0], boundary="bed")
nFaces = len(Xfaces)

Xedges = np.zeros(nFaces + 1)
Xedges[1:-1] = 0.5 * (Xfaces[1:] + Xfaces[:-1])
Xedges[0] = Xfaces[0] - (Xedges[1] - Xfaces[0])
Xedges[-1] = Xfaces[-1] + (Xfaces[-1] - Xedges[-2])

dX = Xedges[1:] - Xedges[:-1]

# Masse suspendue initiale
Cs0 = rdf.readscalar("./", timeList[0], "Cs")
V0 = np.sum(Cs0) * 0.0  # juste pour ne rien supposer ici

# On connaît le volume initial du domaine
Vdomain0 = 1e-4

Ms0 = rhoS * 0.05 * Vdomain0

# Fond initial
zb0 = rdf.readmesh(
    "./", timeList[0], boundary="bed"
)[2]

# Fond final
zbf = rdf.readmesh(
    "./", timeList[-1], boundary="bed"
)[2]

# Masse déposée estimée par le déplacement du fond
dzb = zbf - zb0

Mb = rhoS * CsMax * np.sum(dzb * dX * dy)

print("======================================")
print("Initial suspended mass =", Ms0, "kg")
print("Final deposited mass   =", Mb, "kg")
print("Ratio Mb / Ms0         =", Mb / Ms0)
print("Mass error              =", (Mb - Ms0) / Ms0 * 100, "%")
print("======================================")
