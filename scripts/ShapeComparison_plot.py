################################################################################
 # This file is part of CosTuuM
 # Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 #
 # CosTuuM is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option) any
 # later version.
 #
 # CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 # A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 # details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ###############################################################################

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import CosTuuM

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 5.0)

theta = np.linspace(0.0, 0.5 * np.pi, 50)
ds = [0.25, 0.5, 1.4, np.sqrt(2.0), 1.6, 2.0, 3.0, 4.0, 6.0, -1.0]
dlabels = ["1/4", "1/2", "1.4", "\sqrt{2}", "1.6", "2", "3", "4", "6", "1"]
yax0offsets = [0, 0.3, -4.0, -1.5, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0]
yax1offsets = [-0.02, -0.02, -0.07, -0.02, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0]
colors = pl.cm.ScalarMappable(cmap="tab20").to_rgba(np.linspace(0.0, 1.0, 20))

fig, ax = pl.subplots(2, 1, sharex=True)

results = np.fromfile("ShapeComparison.dat").reshape((len(theta), len(ds), 2))
for id in range(len(ds)):
    d = ds[id]

    quad = results[:, id, :]

    label = None
    if d == -1:
        label = "CDE2"
        color = "k"
    else:
        color = colors[id]
    ax[0].semilogy(theta, quad[:, 0], label=label, color=color)
    ax[1].plot(theta, np.abs(quad[:, 1]) / quad[:, 0], label=label, color=color)
    if d > 0:
        ax[0].text(
            theta[-1] + 0.05,
            quad[-1, 0] + yax0offsets[id],
            "$d={0}$".format(dlabels[id]),
            color=color,
            fontsize=9,
        )
        ax[1].text(
            theta[-1] + 0.05,
            np.abs(quad[-1, 1]) / quad[-1, 0] + yax1offsets[id],
            "$d={0}$".format(dlabels[id]),
            color=color,
            fontsize=9,
        )

ax[0].set_ylabel("$\\langle{}K_{a,I}\\rangle{}$ (code units)")
ax[1].set_ylabel("$P_L$")
ax[1].set_xlabel("$\\theta{}$")
ax[1].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[1].set_xlim(0.0, 1.9)
ax[0].set_ylim(24.0, None)
ax[1].set_ylim(None, 0.75)
ax[1].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])

ax[0].legend(loc="best")
pl.tight_layout()
pl.savefig("ShapeComparison.png", dpi=300, bbox_inches="tight")
pl.savefig("ShapeComparison.pdf", dpi=300, bbox_inches="tight")
