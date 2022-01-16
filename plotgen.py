import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d

from AICC import computeAICC


def getP(h2, h2o):
    r = {}
    r["H2"] = h2
    r["H2O"] = h2o
    r["O2"] = (1 - r["H2"] - r["H2O"]) / 4.76
    if r["O2"] < 0:
        r["O2"] = 0
    if r["H2"] + r["H2O"] > 1:
        r["H2"] = 0
    r["N2"] = 3.76 * r["O2"]
    r["T"] = 300
    r["P"] = 1.01  # in bar
    r["V"] = 1
    res = computeAICC(r)
    if res["P"] > 10000:
        pass
        # print(res)
        # print(r)
    return res["P"]


x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)

X, Y = np.meshgrid(x, y)
Z = np.zeros(shape=(x.shape[0], y.shape[0]))

for i in range(x.shape[0]):
    for j in range(y.shape[0]):
        Z[i, j] = getP(X[i, j], Y[i, j])


fig = plt.figure()

# syntax for 3-D plotting
ax = plt.axes(projection="3d")
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="plasma", edgecolor="none")
ax.set_xlabel("H2")
ax.set_ylabel("H2O")
ax.set_zlabel("Pressure")

plt.show()
