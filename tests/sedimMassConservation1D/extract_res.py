"""

"""

from os import popen
import os.path
import numpy as np
from fluidfoam import readof as rdf


prec = 8

filePath = "./dataSed.txt"


def r(x):
    return str(round(x, prec))


Xmesh, Ymesh, Zmesh = rdf.readmesh("./")
ncells = len(Xmesh)

foamTimes = popen("foamListTimes -withZero").read()
timeList = foamTimes.split("\n")[:-1]
ntimes = len(timeList)

CsFields = np.zeros((ntimes, ncells))

for i, t in enumerate(timeList):
    cs = rdf.readscalar(
        "./", time_name=t, name="Cs", verbose=False)
    if cs.shape == (1,):
        cs = cs * np.ones_like(Xmesh)
    CsFields[i, :] = cs[:]
CsFields = np.round(CsFields, prec)

if os.path.isfile(filePath):
    print(
        f"\nfile {filePath} already exists, "
        + "remove it before saving data")
else:
    with open(filePath, "w") as data:
        data.write("# z; cs(t=" + ") ;cs(t=".join(timeList) + ")\n")
        for i, z in enumerate(Zmesh):
            CsValues = ";".join([str(cs) for cs in CsFields[:, i]])
            line = r(z) + ";" + CsValues + "\n"
            data.write(line)
