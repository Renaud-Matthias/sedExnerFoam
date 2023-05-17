"""
Compare sedimentation simulation results with data stored
in data.txt files

Return error if results and data do not match
"""

import numpy as np
from fluidfoam import readof as rdf

Ysimu = rdf.readmesh('./')[1]

timeList = ['0', '500', '1000', '1500']

Csimu = []

for t in timeList:
    Csimu.append(rdf.readscalar('./', time_name=t, name='C'))

Ydata, c0, c500, c1000, c1500 = np.loadtxt(
    'data.txt', delimiter=';', unpack=True)
Cdata = [c0, c500, c1000, c1500]


tol = 1e-6
success = True

for cs, cd in zip(Csimu, Cdata):
    err = np.max(np.abs(cs - cd) / np.std(cd))
    if err > tol:
        success = False

assert success
