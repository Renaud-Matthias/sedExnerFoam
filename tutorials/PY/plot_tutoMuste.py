"""
Compare experiment results from Muste & Yu (2005) with numerical results
obtained with suspensionFoam solver
Three different cases can be chosen :
- CW : a clear water case with no sediment
- NS1 : a case with natural sand with diameter ranging from 0.21 to 0.25 mm
- NSB1 : a case with crushed nylon particles with very low falling velocity
"""

import os as os
import numpy as np
from scipy import optimize
from scipy import signal

import fluidfoam

from pylab import *
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
from matplotlib import rc
import matplotlib as mpl

rc("text", usetex=True)

# choose case, comment and uncomment

# case = 'CW'
case = "NS1"
# case = 'NBS1'


# Rouse profile
def RouseProfile(alphaRef, z, zRef, Ro):
    alphaRouse = alphaRef * ((1 - z) / z * zRef / (1 - zRef)) ** (Ro)
    return alphaRouse


# Muste data reading
H_f = 0.021
u_star_Muste = 0.042
nu = 1e-6
rhof = 1e3


# *** case selection *** #

solver = "suspensionFoam"

iplotMeanProf = 1
alphaLog = 1
iplotWallUnitProf = 1

if case == "CW":
    sol = "../suspensionFoam/RAS/MUSTE/MUSTE_CW"
    tout = "200/"
    figName = "MUSTE_CW"
    #
    Vel_Muste, Pos_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_Uf.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_UW, UW_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_uw.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_urms, U_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_URMS.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_wrms, W_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_WRMS.txt", unpack=True, delimiter=",", skiprows=0
    )

    alphaMin = 1e-6
    alphaMax = 1e-1
    Wfall = 0.024
    alphaRef = 0.0
    beta = 1
    legloc = 2
elif case == "NS1":
    sol = "../suspensionFoam/RAS/MUSTE/MUSTE_NS1_topoSet"
    tout = "200/"
    figName = "MUSTE_NS1"

    Vel_Muste, Pos_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_NS1_Uf.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_UW, UW_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_uw.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_urms, U_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_URMS.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_wrms, W_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_WRMS.txt", unpack=True, delimiter=",", skiprows=0
    )
    Phi_Muste_Phi, Pos_Muste_Phi = np.loadtxt(
        "../DATA/MUSTE/Phi_NS1.txt", unpack=True, delimiter=",", skiprows=0
    )

    alphaMin = 1e-5
    alphaMax = 2e-2
    legloc = 3
    Wfall = 0.024
    alphaRef = 2e-3
    beta = 2.11
elif case == "NBS1":
    sol = "../suspensionFoam/RAS/MUSTE/MUSTE_NBS1"
    tout = "200/"
    figName = "MUSTE_NBS1"

    Vel_Muste, Pos_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_Uf.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_UW, UW_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_uw.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_urms, U_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_URMS.txt", unpack=True, delimiter=",", skiprows=0
    )
    Pos_Muste_wrms, W_rms_Muste = np.loadtxt(
        "../DATA/MUSTE/CW1_WRMS.txt", unpack=True, delimiter=",", skiprows=0
    )
    Phi_Muste_Phi, Pos_Muste_Phi = np.loadtxt(
        "../DATA/MUSTE/Phi_NB1.txt", unpack=True, delimiter=",", skiprows=0
    )

    alphaMin = 9e-5
    alphaMax = 1.1e-3
    legloc = 2
    Wfall = 0.0006
    alphaRef = 5.2e-4
    beta = 0.52


W_rms_MusteI = np.interp(Pos_Muste_urms, Pos_Muste_wrms, W_rms_Muste)
TKE_Muste = 0.5 * (U_rms_Muste**2 + W_rms_MusteI**2)

dp = 230e-6

# read cordinate matrix
x, z, y = fluidfoam.readmesh(sol)
Z_vect = z
N = np.size(z)

# Read Fields
if case != "CW":
    alpha = fluidfoam.readscalar(sol, tout, "C")
else:
    alpha = np.zeros(N)
k = fluidfoam.readscalar(sol, tout, "k")
Vel = fluidfoam.readvector(sol, tout, "U")
Tauf = fluidfoam.readtensor(sol, tout, "Tau")

u_star_num = (Tauf[3, 0]) ** 0.5

print("u*num=", u_star_num, " u*Muste=", u_star_Muste)
print(r"$z^+=$", Z_vect[0] * u_star_num / nu)
Ro = Wfall / (0.41 * u_star_Muste)
print("Ro=", Ro, " Ro (fit)=", Ro / beta)
#
# linewidth
#
lw = 3
ms = 8
matplotlib.rcParams.update({"font.size": 20})
mpl.rcParams["lines.linewidth"] = lw
mpl.rcParams["lines.markersize"] = ms

#
gs = gridspec.GridSpec(1, 1)
gs.update(
    left=0.065, right=0.975, top=0.95, bottom=0.15, wspace=0.3, hspace=0.3)
gs3 = gridspec.GridSpec(1, 3)
gs.update(
    left=0.065, right=0.975, top=0.95, bottom=0.15, wspace=0.3, hspace=0.3)
gs4 = gridspec.GridSpec(1, 4)
gs.update(
    left=0.065, right=0.975, top=0.95, bottom=0.15, wspace=0.05, hspace=0.15)
gs2 = gridspec.GridSpec(1, 2)
gs.update(
    left=0.065, right=0.975, top=0.95, bottom=0.15, wspace=0.3, hspace=0.3)
#
# Figure size
#
figwidth = 16
figheight = 6

zmin = min(Z_vect / H_f)
zmax = max(Z_vect / H_f)

ExpKeep = np.where(Vel_Muste / u_star_Muste > 11.0)
zRouse = np.linspace(0.05, 0.95, 100)

if iplotMeanProf == 1:
    #
    # plot in dimensioned variables
    #
    fig1 = figure(
        num=1, figsize=(figwidth, figheight),
        dpi=60, facecolor="w", edgecolor="w")

    ax0 = subplot(gs4[0, 0])
    p00 = ax0.plot(Vel[0, :], Z_vect / H_f, "-", label=solver)
    p00 = ax0.plot(
        Vel_Muste[ExpKeep], Pos_Muste[ExpKeep], "or",
        fillstyle="none", label="CW1")
    xlabel(r"$U (m/s)$")
    ylabel(r"$z/h$")
    axis([0, 1, zmin, zmax])

    ax0 = subplot(gs4[0, 1])
    p00 = ax0.plot(
        (1 - alpha[1:N]) * Tauf[3, 1:N],
        Z_vect[1:N] / H_f, "-", label="suspensionFoam")
    p00 = ax0.plot(
        UW_Muste * u_star_Muste**2,
        Pos_Muste_UW, "or", fillstyle="none", label="CW1")
    xlabel(r"$-\rho^f \overline{u^{f\prime} w^{f\prime}}$")
    axis([0, 1.2 * u_star_Muste**2, zmin, zmax])
    ax0.tick_params(labelleft="off")

    ax0 = subplot(gs4[0, 2])
    p00 = ax0.plot(k[:], Z_vect / H_f, "-", label=solver)
    p00 = ax0.plot(
        TKE_Muste * u_star_Muste**2,
        Pos_Muste_urms * nu / (u_star_Muste * H_f),
        "or",
        fillstyle="none",
        label="CW1",
    )
    xlabel(r"$TKE$")
    axis([0, 7 * u_star_Muste**2, zmin, zmax])
    ax0.tick_params(labelleft="off")

    ax0 = subplot(gs4[0, 3])
    p00 = ax0.plot(alpha[:], Z_vect / H_f, "-", label="suspensionFoam")
    p00 = ax0.plot(
        RouseProfile(alphaRef, zRouse, 0.05, Ro), zRouse, "--", label=r"Rouse"
    )
    p00 = ax0.plot(
        RouseProfile(alphaRef, zRouse, 0.05, Ro / beta),
        zRouse,
        "-.",
        label="Rouse (fit)",
    )
    if case == "NS1":
        p00 = ax0.plot(
            Phi_Muste_Phi, Pos_Muste_Phi, "o",
            fillstyle="none", label="NS1")
    if case == "NBS1":
        p00 = ax0.plot(
            Phi_Muste_Phi, Pos_Muste_Phi, "o",
            fillstyle="none", label="NBS1")
    legend(fontsize=12, loc=legloc)
    xlabel(r"$\phi$")
    axis([alphaMin, alphaMax, zmin, zmax])
    if alphaLog == 1:
        ax0.set_xscale("log")
    ax0.tick_params(labelleft="off")

    savefig(
        "Figures/" + figName + "ProfilesDim" + ".eps",
        facecolor="w",
        edgecolor="w",
        format="eps",
    )
    print("save fig")

if iplotWallUnitProf == 1:
    fig2 = figure(num=2, figsize=(12, 8), dpi=60, facecolor="w", edgecolor="w")

    zplus = np.log(np.size(Z_vect))
    ulog = np.zeros(np.size(Z_vect))
    zplus = Z_vect * u_star_Muste / 1e-6
    ulog = 1.0 / 0.41 * np.log(zplus) + 5.5
    toto = np.where(zplus < 11.6)
    ulog[toto] = zplus[toto]

    semilogx(
        Z_vect * u_star_num / nu, Vel[0, :] / u_star_num, "-", label=solver)
    semilogx(
        np.array(Pos_Muste[ExpKeep]) * H_f * u_star_Muste / 1e-6,
        np.array(Vel_Muste[ExpKeep]) / u_star_Muste, "or",
        fillstyle="none", label="Exp")
    semilogx(zplus, ulog, "--", label="Log law")
    ylabel(r"$u/u*$")
    xlabel(r"$z+$")
    axis([1, 1e3, 0, 25])
    legend()


show()
