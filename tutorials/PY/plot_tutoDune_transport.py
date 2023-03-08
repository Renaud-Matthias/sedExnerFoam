"""

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from fluidfoam import readof as rdf

path_tuto = "../exnerFoam/duneTransport"

Z0dune = 0.2
x0dune = 0.3
widthDune = 0.1

timeList = [0, 0.5, 1]  # times at which solution is wanted

# convert timeList to list of string
for i, t in enumerate(timeList):
    if not isinstance(t, str):
        timeList[i] = str(t)

# parameters for Qb, must correspond to exnerFoam case
alpha = 0.05
beta = 2

Qflow = 1.0  # discharge
Hwater = 1.0  # water depth

# ---- load simulation results -----#

X = rdf.readmesh(path_tuto)[0]

ZB_num = []

for t in timeList:
    zb = rdf.readfield(path_tuto, time_name=t, name="Zbvf", boundary="bottom")
    ZB_num.append(zb)

# ---- semi-analytical solution ----#


def CI_dune(X):
    """Initial dune"""
    return Z0dune * np.exp(-(((X - x0dune) / widthDune) ** 2))


def Qb(Zb):
    """bedload flux"""
    return alpha * (Qflow / (Hwater - Zb)) ** beta


def C(Zb):
    """caracteristic equation, dQb/dZb"""
    return (alpha * beta * Qflow**beta) / (Hwater - Zb) ** (beta + 1)


ZB0 = CI_dune(X)

ZB_ana = []

for t in timeList:
    t_float = float(t)
    print(f"\nsolving at t = {t} s")

    def fcost(ZB):
        return np.sum((ZB - CI_dune(X - t_float * C(ZB))) ** 2)

    res = minimize(fcost, ZB0, method="SLSQP", options={"maxiter": 500})
    if res.success:
        print(f"minimization algorithm successfull in {res.nit} iterations")
    else:
        print("minimize algorithm failed after {res.nit} iterations")
        print(f"{res.message}")
    print("precision achieved :")
    print(f"objective function : {res.fun}")
    print(f"jacobian value : {res.fun}")
    ZB_ana.append(res.x)


def RMSE_SCORE(Xref, X):
    """return RMSE and Score based RMSE"""
    rmse = np.sqrt(np.mean((Xref - X) ** 2))
    score = 1 - rmse / np.std(Xref)
    return rmse, score


# error between analytical and numerical solutions
err_rmse = []
err_score = []
for zb_ana, zb_num in zip(ZB_ana, ZB_num):
    rmse, score = RMSE_SCORE(zb_ana, zb_num)
    err_rmse.append(rmse)
    err_score.append(score)

# PLOT analytical solution #

fig1, axs = plt.subplots(3, figsize=(10, 10), layout="constrained")

i = 0
for zb_ana, zb_num in zip(ZB_ana, ZB_num):
    axs[0].plot(X, zb_ana, label=f"t = {timeList[i]} s")
    axs[1].plot(X, zb_num)
    axs[2].plot(X, np.abs(zb_ana - zb_num))
    i += 1
# set label
axs[0].set_title("semi-analytical solution")
axs[1].set_title("numerical solution")
axs[2].set_title("error absolute")

axs[0].set_xlabel("x")
axs[0].set_ylabel("zb")
axs[1].set_xlabel("x")
axs[1].set_ylabel("zb")
axs[2].set_xlabel("x")
axs[2].set_ylabel(r"$|zb_{ref} - zb_{num}|$")

axs[0].legend()

if len(timeList) > 1:
    fig2, ax = plt.subplots()
    t_err = [float(t) for t in timeList]  # corresponding time
    ax.plot(t_err, err_score, marker="o")
    ax.axhline(y=1, label="perfect score", color="black")

    ax.legend()

plt.show()
