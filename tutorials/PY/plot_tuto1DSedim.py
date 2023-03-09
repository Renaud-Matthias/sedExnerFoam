import numpy as np
from fluidfoam import readof as rdf
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

path_data = "../DATA/PhamVanBang2006/"
path_tuto = "../suspensionFoam/laminar/1DSedim/"

t_plot_list = [180, 540, 1140, 1500]
# possible choices are :
# 0.0 ; 60.0 ; 120.0 ; 180.0 ; 240.0 ; 300.0 ; 360.0 ; 420.0 ; 480.0 ; 540.0
# 600.0 ; 660.0 ; 720.0 ; 780.0 ; 840.0 ; 900.0 ; 960.0 ; 1020.0 ; 1080.0
# 1140.0 ; 1200.0 ; 1260.0 ; 1320.0 ; 1380.0 ; 1440.0 ; 1500.0 ; 1560.0
# 1620.0 ; 1680.0 ; 1740.0

# --- READ NUMERICAL RESULTS --- #

Ynum = rdf.readmesh(path_tuto)[1] - 0.1

Cnum = []
for t in t_plot_list:
    Cnum.append(rdf.readscalar(path_tuto, time_name=str(t), name="C"))

# --- READ EXPERIMENTAL RESULTS --- #


def read_exp(path_file):
    time_list = []
    z_dict = {}
    c_dict = {}
    with open(path_file, "r") as f:
        for line in f:
            if "#" in line:
                continue
            t, z, c = line.split(";")
            t = float(t)
            z = np.array([float(ch) for ch in z.split(" ")])
            c = np.array([float(ch) for ch in c.split(" ")])
            time_list.append(t)
            z_dict.update({t: z})
            c_dict.update({t: c})
    return (time_list, z_dict, c_dict)


# load experimental results, concentration profiles
time_list, z_dict, c_dict = read_exp(path_data + "concentration.txt")
print("times at which experimental data are available (in seconds) :")
print(" ; ".join([str(t) for t in time_list]))
z_array = z_dict[0]  # possible because all z array are identical

# load experimental results, interface position
t_int, y_int_sup, y_int_inf = np.loadtxt(
    path_data + "./interfaces.txt", unpack=True)

# --- CORRECT TIME EXPERIMENT --- #


t_end = 1040  # time at which interfaces join

t_fit = []  # time liste for linear fit
for t in t_int:
    if t < t_end:
        t_fit.append(t)
t_fit = np.array(t_fit)
n_elements = len(t_fit)
y_new_inf = y_int_inf[:n_elements]

# linear regression coefficients
slope_fit = (np.mean(t_fit * y_new_inf)
             - np.mean(t_fit)
             * np.mean(y_new_inf)) / np.var(t_fit)

y0_fit = np.mean(y_new_inf) - slope_fit * np.mean(t_fit)

t0_true = (-0.1 - y0_fit)/slope_fit
print("real start time of experiment : ", round(t0_true, 3), " s")

# correct t_int
t_int = t_int - t0_true


def initial_Cexp(Z):
    """initial volumic concentration profile"""
    C_profile = 0.5 * 0.5 * (1 + np.tanh((Z - 0.054) / (0.049 - Z)*10))
    return np.where(Z > 0.049, C_profile, 0.5)


def interp_Cexp(t_val):
    """Get value of c at time t_val through interpolation"""
    if t_val < 0:
        t1, t2 = t0_true, time_list[0]
        c1, c2 = initial_Cexp(z_dict[t2] + 0.1), c_dict[t2]
    else:
        for i in range(len(time_list)):
            if time_list[i+1] > t_val:
                break
        t1, t2 = time_list[i], time_list[i+1]
        c1, c2 = c_dict[t1], c_dict[t2]
    dt1, dt2 = t_val - t1, t2 - t_val
    dt = dt1 + dt2
    c_val = ((dt - dt1) * c1 + (dt - dt2) * c2) / dt
    return c_val


# --- FIND INTERFACE POSITION --- #


def find_pos_inter_inf(Carr, Yarr, cmax=0.6, tol=0.05):
    """Return position of lower interface"""
    c_inter = cmax * (1 - tol)
    for c, y in zip(Carr, Yarr):
        if c < c_inter:
            return y
    print("no interface found")
    return None


def find_pos_inter_sup(Carr, Yarr, cmax=0.6, tol=0.05):
    """Return position of upper interface"""
    c_inter = cmax * tol
    for c, y in zip(Carr, Yarr):
        if c < c_inter:
            return y
    print("no interface found")
    return None


Yinter_inf = np.zeros(len(time_list))
Yinter_sup = np.zeros(len(time_list))

for i, t in enumerate(time_list):
    Ctemp = rdf.readscalar(path_tuto, time_name=str(round(t)), name="C")
    Yinter_inf[i] = find_pos_inter_inf(Ctemp, Ynum)
    Yinter_sup[i] = find_pos_inter_sup(Ctemp, Ynum)

# --- PLOT RESULTS --- #


nplot = len(Cnum)
nlines = 1 + (nplot - 1) // 2
ncolumns = min(2, nplot)

fig1 = plt.figure(figsize=(4 * ncolumns, 4 * nlines), layout="constrained")
col_exp, col_num = "black", "#ea5545"

gs = fig1.add_gridspec(nlines, ncolumns)
axs = []

for i, t in enumerate(t_plot_list):
    i_line, i_col = divmod(i, 2)
    axs.append(fig1.add_subplot(gs[i_line, i_col]))
    axs[i].plot(
        interp_Cexp(t + t0_true), z_array * 100 + 10,
        ls="--", lw=2, c=col_exp, label="experiment")
    axs[i].plot(
        Cnum[i], Ynum * 100 +10, lw=2, c=col_num, label="OpenFoam")

    axs[i].set_title(f"t = {t} s")
    axs[i].set_ylim(0, 6)
    # add y label only if plot is on left side
    if i_col == 0:
        axs[i].set_ylabel("y [cm]")
    if i_line == nlines - 1:
        axs[i].set_xlabel("c")
        if i_col == 0:
            axs[i].legend(loc="lower left")

fig2 = plt.figure(figsize=(6, 4))
col_up, col_lo = "#1984c5", "#a57c1b"

ax = fig2.add_subplot()
ax.plot(t_int, y_int_sup * 100 + 10,
        marker="o", ls="None", c=col_up)
ax.plot(t_int, y_int_inf * 100 + 10,
        marker="o", ls="None", c=col_lo)
ax.plot(time_list, Yinter_sup * 100 + 10,
        c=col_up, lw=2)
ax.plot(time_list, Yinter_inf * 100 + 10,
        c=col_lo, lw=2)

ax.set_xlabel("t [s]")
ax.set_ylabel("z [cm]")

legend_elements = [
    Patch(facecolor=col_up, label="upper interface"),
    Patch(facecolor=col_lo, label="lower interface"),
    Line2D([0], [0], marker="o", ls="None",
           color="black", label="experimental"),
    Line2D([0], [0], color="black", label="numerical")
]

ax.legend(
    handles=legend_elements, loc="lower right",
    prop={"size": 10}, frameon=False)

plt.show()
