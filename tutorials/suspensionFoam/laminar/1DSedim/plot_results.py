"""
"""

import matplotlib.pyplot as plt
import numpy as np
from fluidfoam import readof as rdf

nu_f = 1e-6  # fluid viscosity

X, Y, Z = rdf.readmesh("./")

# timeList = [1,2,3,4,5]#,180,240,480]
timeList = [180, 240, 480]
timeStr = [str(tvalue) for tvalue in timeList]
nt = len(timeStr)

C_field = []
for t in timeStr:
    C_field.append(rdf.readscalar("./", time_name=t, name="C"))

Ws_field = []
for t in timeStr:
    Ws_field.append(rdf.readfield("./", time_name=t, name="Ws")[1])


fig = plt.figure(figsize=(12, 6), layout="constrained")
spec = fig.add_gridspec(1, 2)

ax0 = fig.add_subplot(spec[0, 0])
ax1 = fig.add_subplot(spec[0, 1])

for C, i in zip(C_field, range(nt)):
    ax0.plot(C, Y, label=f"t = {timeStr[i]} s")
ax0.set_ylabel("y [m]")
ax0.set_xlabel("C")
ax0.set_xlim(None, 1)
ax0.legend()

for Ws, i in zip(Ws_field, range(nt)):
    ax1.plot(-Ws, Y, label=f"t = {timeStr[i]} s")
ax1.set_ylabel("y [m]")
ax1.set_xlabel("Ws [m/s]")

plt.show()
