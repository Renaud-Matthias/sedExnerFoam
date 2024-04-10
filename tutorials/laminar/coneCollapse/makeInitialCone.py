"""
create initial dune geometry
"""

from math import exp, pi, tan
import os

path = "./constant/polyMesh"

# domain dimensions
domHeight = 0.5

angleCone = 40.  # degrees

coneRadius = 0.5
coneHeight = coneRadius * tan(angleCone * pi / 180)

if coneHeight > domHeight:
    print("cone height is larger than domain height")
    raise ValueError("domain height must be larger than cone height")


def zCone(x):
    zout = 0.
    if x < coneRadius:
        zout = coneHeight * (1 - (1/coneRadius) * x)
    return zout


def newZcoord(x, z):
    zc = zCone(x)
    zout = domHeight**2 - (domHeight - z) * (domHeight - zc)
    return zout/domHeight


newlines = []

os.system(f"cp {path}/points {path}/newPoints")

newPoints = open(f"{path}/newPoints", "w")

with open(f"{path}/points", "r") as f:
    for line in f:
        newline = line
        if line[0] == "(" and line[-2] == ")":
            xs, ys, zs = line[:-1].strip("()").split(" ")
            # print(xs, ys, zs)
            x, z = float(xs), float(zs)
            z = newZcoord(x, z)
            zs = str(round(z, 10))
            newline = f"({xs} {ys} {zs})\n"
        newPoints.write(newline)

newPoints.close()

os.system(f"rm {path}/points && mv {path}/newPoints {path}/points")
