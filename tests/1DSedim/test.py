"""
Compare sedimentation simulation results with data stored
in data.txt files

Return error if results and data do not match
the data is the results of a previous simulation
"""

import numpy as np
from fluidfoam import readof as rdf
import os

print(" --- running test 1DSedimentation (MULES formulation) --- ")

success = True
tol = 1e-6

Ysimu = rdf.readmesh('./', verbose=False)[1]

foamTimes = os.popen("foamListTimes").read()
timeList = foamTimes.split("\n")[:-1]

# load results from previous simulation
Ydata, *csData = np.loadtxt(
    "dataSedim.txt", delimiter=";", unpack=True)

for i, t in enumerate(timeList):
    Cs = rdf.readscalar(
        "./", time_name=t, name="Cs", verbose=False)
    
    err = (Cs - csData[i])
    if np.any(err > tol):
        success = False
        print(
            "ERROR! results not matching previous simulation\n"
            + f"maximum error on volume fraction: {np.max(err)}")

assert success
