"""
Compare experiment results from Muste & Yu (2005) with numerical results
obtained with suspensionFoam solver
Three different cases can be chosen :
- CW : a clear water case with no sediment
- NS1 : a case with natural sand with diameter ranging from 0.21 to 0.25 mm
- NSB1 : a case with crushed nylon particles with very low falling velocity
Tutorial rouse profile with validation thanks to Muster experiment
"""

import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readof as rdf

pathFoam = "../RAS/MUSTE/"

# plot parameters
col_exp = "#F05039"
col_num = "#3D65A5"
col_num_bis = "#2DB085"
col_rouse = "black"
markersize = 5

# Muste data
H_water = 0.021  # water height
u_star_muste = 0.042  # experimental friction velocity
nu = 1.e-6  # water kinematic viscosity
rhof = 1e3  # water density
WfallNS1 = 0.024  # fall velocity of natural sand
WfallNBS1 = 6e-4  # fall velocity of crushed nylon particle


def RouseProfile(Cref, z, zRef, Ro):
    """Returns the Rouse concentration profile along z"""
    C_Rouse = Cref * (((H_water - z) / z) * (zRef / (H_water - zRef))) ** (Ro)
    return C_Rouse


# NS1 case #

print("Reading NS1 case")
C_MusteNS1, Z_MusteNS1 = np.loadtxt(
        "../DATA/MUSTE/Phi_NS1.txt", unpack=True, delimiter=",", skiprows=0)

betaNS1 = 2.11  # sediment and momentum diffusion coefficient
RoNS1 = WfallNS1 / (0.41 * betaNS1 * u_star_muste)  # Rouse number
idxCmaxNS1 = np.argmax(C_MusteNS1)
crefNS1 = C_MusteNS1[idxCmaxNS1]  # reference concentration
zrefNS1 = Z_MusteNS1[idxCmaxNS1] * H_water  # reference height
Zrouse = np.linspace(1e-5, 0.95 * H_water, 100)
CrouseNS1 = RouseProfile(crefNS1, Zrouse, zrefNS1, RoNS1)
print("Rouse number : ", RoNS1)
print(crefNS1, zrefNS1)

OFNS1case = pathFoam + "MUSTE_NS1"
Z_foamNS1 = rdf.readmesh(OFNS1case)[1]
C_foamNS1 = rdf.readscalar(OFNS1case, "latestTime", "Cs")

# same as NS1 but with reference concentration at bottom
OFNS1TScase = pathFoam + "MUSTE_NS1_topoSet"
Z_foamNS1TS = rdf.readmesh(OFNS1TScase)[1]
C_foamNS1TS = rdf.readscalar(OFNS1TScase, "latestTime", "Cs")

cminNS1 = 1e-5
cmaxNS1 = 5e-2

# NBS1 case #

print("\nReading NBS1 case")
C_MusteNBS1, Z_MusteNBS1 = np.loadtxt(
        "../DATA/MUSTE/Phi_NB1.txt", unpack=True, delimiter=",", skiprows=0)

betaNBS1 = 0.52
RoNBS1 = WfallNBS1 / (0.41 * betaNBS1 * u_star_muste)
idxCmaxNBS1 = np.argmax(C_MusteNBS1)
crefNBS1 = C_MusteNS1[idxCmaxNBS1]
zrefNBS1 = Z_MusteNBS1[idxCmaxNBS1] * H_water
crefNBS1 = max(C_MusteNBS1)
zrefNBS1 = 0.05 * H_water
Zrouse = np.linspace(1e-5, 0.95 * H_water, 100)
CrouseNBS1 = RouseProfile(crefNBS1, Zrouse, zrefNBS1, RoNBS1)
print("Rouse number : ", RoNBS1)

OFNBS1case = pathFoam + "MUSTE_NBS1"
Z_foamNBS1 = rdf.readmesh(OFNBS1case)[1]
C_foamNBS1 = rdf.readscalar(OFNBS1case, "latestTime", "Cs")

cminNBS1 = 1e-4
cmaxNBS1 = 2e-3


# plot

fig = plt.figure(figsize=(10, 5), layout="constrained")
gs = fig.add_gridspec(1, 2)
axNS1 = fig.add_subplot(gs[0, 0])

axNS1.plot(
    C_MusteNS1, Z_MusteNS1, marker="o", markersize=markersize,
    ls="none", markerfacecolor=col_exp, markeredgecolor="black",
    markeredgewidth=0.5, label="experiment")
axNS1.plot(
    CrouseNS1, Zrouse / H_water, color=col_rouse,
    ls="dashed", label="Rouse profile")
axNS1.plot(
    C_foamNS1, Z_foamNS1 / H_water, color=col_num,
    lw=2, label="OF, zeroGradient BC")
axNS1.plot(
    C_foamNS1TS, Z_foamNS1TS / H_water, color=col_num_bis,
    lw=2, label="OF, reference level")

axNS1.set_title("NS1 : natural sand")
axNS1.set_xlabel("Concentration")
axNS1.set_ylabel("z/H")
axNS1.set_xscale("log")
axNS1.set_xlim(cminNS1, cmaxNS1)
axNS1.set_ylim(0, 1)
axNS1.grid(True)
axNS1.legend()

axNBS1 = fig.add_subplot(gs[0, 1])

axNBS1.plot(
    C_MusteNBS1, Z_MusteNBS1, marker="o",
    markersize=markersize, ls="none", markerfacecolor=col_exp,
    markeredgecolor="black", markeredgewidth=0.5)
axNBS1.plot(
    CrouseNBS1, Zrouse / H_water, color=col_rouse,
    ls="dashed", label="Rouse profile")
axNBS1.plot(
    C_foamNBS1, Z_foamNBS1 / H_water, color=col_num,
    lw=2, label="OpenFoam")

axNBS1.set_title("NBS1 : crushed nylon")
axNBS1.set_xlabel("Concentration")
axNBS1.set_xscale("log")
axNBS1.set_xlim(cminNBS1, cmaxNBS1)
axNBS1.set_ylim(0, 1)
axNBS1.set_yticklabels([])
axNBS1.grid(True)

plt.show()

answer = ""

while answer.lower() not in ("yes", "no"):
    answer = input("do you want to save the figure? (yes / no)")

if answer.lower() == "yes":
    fig.savefig(
        "./Figures/Muste_concentration_profiles.png",
        format="png", transparent=True)
