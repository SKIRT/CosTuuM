import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True

rad_to_deg = 180.0 / np.pi

# load the data
data = np.loadtxt("test_tmatrixcalculator_result.txt")

# retrieve the grids in theta and phi and the single and ensemble scattering
# coefficients and reshape them to the 100x100 grid in angles
theta = data[:, 0].reshape((100, 100))
phi = data[:, 1].reshape((100, 100))
Zsingle = data[:, 2].reshape((100, 100))
Zensemble = data[:, 3].reshape((100, 100))

# create two axes to plot the two grids
fig, ax = pl.subplots(2, 1)

# plot the single particle grid (perfect alignment)
ax[0].contourf(phi * rad_to_deg, theta * rad_to_deg, np.log10(Zsingle), 500)
ax[0].set_ylabel("$\\theta{}$")
ax[0].set_xlabel("$\\phi{}$")
ax[0].set_title("Perfect alignment")

# plot the ensemble grid (imperfect alignment)
ax[1].contourf(phi * rad_to_deg, theta * rad_to_deg, np.log10(Zensemble), 500)
ax[1].set_ylabel("$\\theta{}$")
ax[1].set_xlabel("$\\phi{}$")
ax[1].set_title("Imperfect alignment")

# save the figure
pl.tight_layout()
pl.savefig("test_tmatrixcalculator_result.png", dpi=300, bbox_inches="tight")
