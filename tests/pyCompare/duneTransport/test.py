"""

"""

import numpy as np
from fluidfoam import readof as rdf


print(" --- running test dune transport --- ")

time = "latestTime"

success = True
tol = 1e-4

verbose = False

# physical parameters

dS = 0.4e-3  # sediment diameter
rhoS = 2500.  # sediment density
rhoF = 1000.  # fluid density
nuF = 1e-6  # fluid kinematic viscosity
g = 9.81  # gravity acceleration
betaRepDeg = 28.  # granular material repose angle in degrees
betaRep = betaRepDeg * np.pi / 180

einNum = np.sqrt((rhoS/rhoF - 1) * g * dS**3)

Dstar = dS * ((rhoS/rhoF - 1)*g / nuF**2)**(1/3)

critShSoulsby = (0.3 / (1 + 1.2*Dstar)) + 0.055*(1 - np.exp(-0.02*Dstar))

if verbose:
    print(f"dimless diameter, Dstar = {Dstar}")
    print(f"Einstein number, Ei = {einNum} m2/s")
    print(f"critical Shields number, {critShSoulsby}")

Xb, Yb, Zb = rdf.readmesh("./", time, boundary="bed", verbose=False)

critShields = rdf.readscalar(
    "./", time, "critShieldsVf", boundary="bed", verbose=False)
if critShields.shape == (1,):
    critShields = np.ones_like(Xb) * critShields
shX, shY, shZ = rdf.readvector(
    "./", time, "shieldsVf", boundary="bed", verbose=False)
magSh = np.sqrt(shX**2 + shZ**2)

# compute slope angle beta from Shields components
beta = np.arctan(shZ / shX)
cosAlpha = np.where(shX*beta > 0, -1, 1)
# limit beta values  to betaRep
betaEff = np.where(np.abs(beta) < betaRep, np.abs(beta), betaRep)
slopeCorr = np.cos(betaEff) - cosAlpha * np.sin(betaEff)/np.tan(betaRep)
critShCorr = critShSoulsby * slopeCorr

qbX, qbY, qbZ = rdf.readvector(
    "./", time, "qbVf", boundary="bed", verbose=False)
magQb = np.sqrt(qbX**2 + qbZ**2)

# bedload from Meyer-Peter & Müller formula
qbMPM = einNum * 8 * np.where(
    magSh > critShields, magSh - critShields, 0)**1.5
# bedload due to avalanche from Vinent et. al (2019) formula
Qav = 5e-3
qbAv = Qav * np.where(
    np.abs(beta) > betaRep,
    (np.tanh(np.tan(np.abs(beta)))-np.tanh(np.tan(betaRep)))
    / (1-np.tanh(np.tan(betaRep))), 0)
# total bedload, Meyer-Pter + avalanche
qbTot = np.where(beta*shX > 0, np.abs(qbAv-qbMPM), qbAv + qbMPM)

# relative error on critical Shields number
# with slope Correction
errCritShields = np.max((critShields - critShCorr)/np.max(critShCorr))
if errCritShields > tol:
    success = False
    print("error on critical Shields values"
          + f"\n relative error is {errCritShields*100} %"
          + f"\n tolerance is {tol*100} %")

# relative error on bedload
errQb = np.max((magQb - qbTot)/np.max(qbTot))
if errQb > tol:
    success = False
    print("error on bedload values"
          + f"\n relative error is {errQb*100} %"
          + f"\n tolerance is {tol*100} %")

assert success
