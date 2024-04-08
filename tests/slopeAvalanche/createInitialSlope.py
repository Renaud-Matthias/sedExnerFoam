"""
create initial dune geometry
"""

import numpy as np
import os

path = "./constant/polyMesh"

Hz = 0.1  # domain height
Lx = 0.1  # domain length
Ls = 0.1  # length of slope area

slopeDeg = 40  # initial in degrees

newlines = []
slopeRad = slopeDeg * np.pi / 180
dH = np.tan(slopeRad) * Ls  # height difference


def getZbed(x):
    z = dH - np.tan(slopeRad) * (x - 0.5*(Lx - Ls))
    if x < 0.5 * (Lx - Ls):
        z = dH
    elif x > 0.5 * (Lx + Ls):
        z = 0
    return z


os.system(f"cp {path}/points {path}/newPoints")

newPoints = open(f"{path}/newPoints", "w")

with open(f"{path}/points", "r") as f:
    for line in f:
        newline = line
        if line[0] == "(" and line[-2:] == ")\n":
            xs, ys, zs = line[:-1].strip("()").split(" ")
            x, z = float(xs), float(zs)
            zb = getZbed(x)
            znew = zb + z * (Hz-zb)/Hz
            zs = str(round(znew, 10))
            newline = f"({xs} {ys} {zs})\n"
        newPoints.write(newline)

newPoints.close()

os.system(f"rm {path}/points && mv {path}/newPoints {path}/points")
