import numpy as np
import matplotlib.pyplot as plt

path_data = '../DATA/PhamVanBang2006/'

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

time_list, z_dict, c_dict = read_exp(path_data+'concentration.txt')
print('times at which experimental data are available (in seconds) :')
print(' ; '.join([str(t) for t in time_list]))

t_int, z_int1, z_int2 = np.loadtxt(path_data+'./interfaces.txt',unpack=True)

fig = plt.figure(figsize=(9,6))
gs = fig.add_gridspec(1,2)

ax0 = fig.add_subplot(gs[0,0])
ax0.plot(c_dict[1200],z_dict[1200],ls='--')

ax1 = fig.add_subplot(gs[0,1])
ax1.plot(t_int,z_int1,marker='o',ls='None')
ax1.plot(t_int,z_int2,marker='o',ls='None')

plt.show()
