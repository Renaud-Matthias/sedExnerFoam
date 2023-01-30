import numpy as np
from fluidfoam import readof as rdf
import matplotlib.pyplot as plt

path_data = '../DATA/PhamVanBang2006/'
path_tuto = '../suspensionFoam/laminar/1DSedim/'

t_plot_list = [180, 600, 900, 1200]
# possible choices :
# 0.0 ; 60.0 ; 120.0 ; 180.0 ; 240.0 ; 300.0 ; 360.0 ; 420.0 ; 480.0 ; 540.0 ; 600.0 ; 660.0 ; 720.0 ; 780.0 ; 840.0 ; 900.0 ; 960.0 ; 1020.0 ; 1080.0 ; 1140.0 ; 1200.0 ; 1260.0 ; 1320.0 ; 1380.0 ; 1440.0 ; 1500.0 ; 1560.0 ; 1620.0 ; 1680.0 ; 1740.0

#### READ NUMERICAL RESULTS #####

Ynum = rdf.readmesh(path_tuto)[1] - 0.1

Cnum = []
for t in t_plot_list :
    Cnum.append(rdf.readscalar(path_tuto,time_name=str(t),name='C'))

### READ EXPERIMENTAL RESULTS ###

def read_exp(path_file) :
    time_list = []
    z_dict = {}
    c_dict = {}
    with open(path_file,'r') as f :
        for line in f :
            if '#' in line :
                continue
            t,z,c = line.split(';')
            t = float(t)
            z = np.array([float(ch) for ch in z.split(' ')])
            c = np.array([float(ch) for ch in c.split(' ')])
            time_list.append(t)
            z_dict.update({t : z})
            c_dict.update({t : c})
    return (time_list, z_dict, c_dict)

# load experimental results, concentration profiles
time_list, z_dict, c_dict = read_exp(path_data+'concentration.txt')
print('times at which experimental data are available (in seconds) :')
print(' ; '.join([str(t) for t in time_list]))

# load experimental results, interface position
t_int, y_int1, y_int2 = np.loadtxt(path_data+'./interfaces.txt',unpack=True)

### FIND INTERFACE POSITION ###

def find_pos_inter(Carr, Yarr, cmax=0.6, tol=0.02) :
    c_inter = cmax*(1-tol)
    for i in range(len(time_list)-1,-1,-1) :
        if Carr[i] > c_inter :
            return Yarr[i]
    print('no interface found')
    return None

#def fin_pos_inter_sup(Carr, Yarr, cmax=0.6, tol=0.02) :
#    c_inter = cmax*tol
#    for i in range(len(time_list)-1,-1,-1) :
#        if Carr[i] 

Yinter = np.zeros(len(time_list))

for i,t in enumerate(time_list) :
    Ctemp = rdf.readscalar(path_tuto,time_name=str(round(t)),name='C')
    Yinter[i] = find_pos_inter(Ctemp,Ynum,tol=0.1)

Yinter[Yinter == None] = -0.1
                  
print(Yinter)

### PLOT RESULTS ###

nplot = len(Cnum)
nlines = nplot//2 + 2*(nplot%2)

fig1 = plt.figure(figsize=(8,4*nlines))
gs = fig1.add_gridspec(2,nlines)
axs = []

for  i,t in enumerate(t_plot_list) :
    axs.append(fig1.add_subplot(gs[i%2,i//2]))
    axs[i].plot(c_dict[t],z_dict[t],ls='--')
    axs[i].plot(Cnum[i], Ynum)
    axs[i].set_title(f't = {t} s')
    #axs[i].hline(Yinter[

fig2 = plt.figure(figsize=(6,4))
ax = fig2.add_subplot()
ax.plot(t_int,y_int1,marker='o',ls='None')
ax.plot(t_int,y_int2,marker='o',ls='None')
ax.plot(time_list, Yinter)

plt.show()
