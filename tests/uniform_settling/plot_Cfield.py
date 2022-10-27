'''
Plot the concentration field at the last time step
'''

from fluidfoam import readof as rdf
import matplotlib.pyplot as plt

X,Y,Z = rdf.readmesh('./') # read mesh in current folder

C_field = rdf.readscalar('./','latestTime',name='C')

print(X.shape)
print(Z.shape)
print(C_field.shape)

fig, ax = plt.subplots()

ax.set_title('concentration of suspended sediment')
ax.tripcolor(X,Z,C_field,cmap='jet',shading='flat')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')

plt.show()
