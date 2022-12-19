
import matplotlib.pyplot as plt
import numpy as np
from fluidfoam import readof as rdf

nu_f = 1e-6 # fluid viscosity

X,Y,Z = rdf.readmesh('./')

timeList = [0,2,4,6,8,10]
timeStr = [str(tvalue) for tvalue in timeList]
nt = len(timeStr)

C_field = []
Ux,Uy,Uz = rdf.readfield('./',time_name='latestTime',name='U')

for t in timeStr :
    C_field.append(rdf.readscalar('./',time_name=t,name='C'))


fig = plt.figure(figsize=(6,6),layout='constrained')
spec = fig.add_gridspec(1,1)

ax0 = fig.add_subplot(spec[0,0])


for C,i in zip(C_field,range(nt)) :
    ax0.plot(C,Y,label=f't = {timeStr[i]} s')
ax0.set_ylabel('y [m]')
ax0.set_xlabel('C')
ax0.set_xlim(0,0.5)
ax0.legend()

plt.show()
