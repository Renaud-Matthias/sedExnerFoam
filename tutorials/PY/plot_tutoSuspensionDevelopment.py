"""
Plot iso line of suspended sediment concentration
compare with pseudo analytical solution
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from fluidfoam import readof as rdf
import sys
import os
from hjelmfeltLenauSolution import suspensionDevSol


pathCase = "../RAS/suspensionDev"

time = "latestTime"

save = True

# physical parameters
Hwater = 0.1  # water depth
critShields = 0.05  # critical Shields number
g = 9.81  # gravity acceleration, m/s2
dS = 0.12e-3  # sediment diameter, m
rhoS = 2650.  # sediment density
nuF = 1e-6  # fluid kinematic viscosity, m2/s
rhoF = 1000.  # fluid density
kappa = 0.41  # von Karman constant
Sc = 1.  # Schmidt number
aRef = 0.05 * Hwater  # reference level, m

# levels for Cs iso lines plot
levels = [0.05, 0.1, 0.2, 0.3, 0.5]


def CrefVanRijn(zref, shields, critShields):
    """
    Compute reference concentration, van Rijn (1984)
    zref: float reference height, m
    shields: float, shields Number
    critShields: float, critical Shields number
    """
    # dimensionless relative shear stress
    Tshear = shields/critShields - 1
    # dimension less sediment diameter
    Dstar = dS * ((((rhoS/rhoF - 1) * g) / nuF)**2)**(1/3)
    cref = 0.015 * (dS/zref) * Tshear**1.5 * Dstar**(-0.3)
    return cref


Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase)
Xbed, Ybed, Zbed = rdf.readmesh(pathCase, boundary="bed")
xmin, xmax = 0, np.max(Xmesh)

# read cell volumes
try:
    Vcells = rdf.readscalar(pathCase, time, "V")
except FileNotFoundError:
    os.system(f"postProcess -case {pathCase} -func writeCellVolumes")
    Vcells = rdf.readscalar(pathCase, time, "V")
# read U field
Ufield = rdf.readvector(pathCase, time, "U")
# compute Umean
Umean = np.sum(np.linalg.norm(Ufield, axis=0)*Vcells)/np.sum(Vcells)
print(f"- mean velocity magnitude: {Umean} m/s\n")

# shields number
Shields = rdf.readvector(pathCase, time, "shieldsVf", boundary="bed")[0]
shields = np.mean(Shields)
print(
    f"- shields: max = {np.max(Shields)}, min = {np.min(Shields)}, "
    + f"mean = {np.mean(Shields)}\n")
Ustar = np.sqrt(Shields * (rhoS/rhoF - 1) * g * dS)
ustar = np.mean(Ustar)
print(f"- uf: max = {np.max(Ustar)}, min = {np.min(Ustar)}, "
      + f"mean = {np.mean(Ustar)}\n")


# get settling velocity
ws = np.abs(rdf.readvector(pathCase, time, "Ws", verbose=False)[2, 0])
print(f"- settling velocity ws = {ws} m/s")
# rouse number
Ro = ws / (kappa * ustar * Sc)
print(f"- Rouse number = {Ro}")

# get reference concentration
cref = CrefVanRijn(aRef, shields, critShields)
print(f"- reference contration = {cref}")
cref = 0.025

# sedExnerFoam suspended load concentration
CsField = rdf.readscalar(pathCase, time, "Cs")
vmin = 1e-7
CsField = np.where(CsField > vmin, CsField, vmin)

# get pseudo analytical solution, Hjelmfelt & Lenau (1970)
suspDevSol = suspensionDevSol(Ro, aRef / Hwater)
Xgrid = np.linspace(xmin/Hwater, xmax*(Sc*kappa*ustar/(Umean*Hwater)), 200)
Zgrid = np.linspace(aRef/Hwater, 1, 50)
Xgrid, Zgrid = np.meshgrid(Xgrid, Zgrid)

# compute solution
Cpseudo = suspDevSol.get_solution(Xgrid, Zgrid)

Xgrid *= Umean / (kappa * ustar * Sc)


fig, axCs = plt.subplots(figsize=(8.3, 6))

CS = axCs.tricontour(
    Xmesh/Hwater, Zmesh/Hwater, CsField/cref,
    levels=levels, linewidths=3., colors=["cornflowerblue"])

axCs.contour(
    Xgrid, Zgrid, Cpseudo, levels=levels,
    linestyles=["dotted"], linewidths=3., colors=["crimson"])

axCs.set_xlabel("x/H", fontsize=15)
axCs.set_ylabel("z/H", fontsize=15)
axCs.tick_params(
    axis='both', which='major', labelsize=15)
axCs.grid()
# text on plot
axCs.text(10, 0.9, r"$R_o=$"+f"{round(Ro, 2)}", fontsize=20)
axCs.text(140, 0.1, "0.5", fontsize=15)
axCs.text(140, 0.24, "0.3", fontsize=15)
axCs.text(140, 0.4, "0.2", fontsize=15)
axCs.text(140, 0.71, "0.1", fontsize=15)
axCs.text(140, 0.88, "0.05", fontsize=15)

fig.tight_layout()

plt.show()

if save:
    imname = "suspensionDev_sedExnerFoam_Ro_0_5.eps"
    fig.savefig("./" + imname, format="eps", transparent=True)
