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

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 3.0)

thetas = np.linspace(0.0, 0.5 * np.pi, 50)

nmin = 1
nmax = 20
nmax_plot = 6
cm = pl.get_cmap("viridis")
cb = pl.contourf(
    [[0.0, 0.0], [0.0, 0.0]], np.linspace(nmin, nmax_plot, 500), cmap=cm
)
pl.clf()

colwidth = 5
ax = [
    pl.subplot2grid((2, colwidth + 1), (0, 0), colspan=colwidth),
    pl.subplot2grid((2, colwidth + 1), (1, 0), colspan=colwidth),
    pl.subplot2grid((2, colwidth + 1), (0, colwidth), rowspan=2),
]
ref_prolate = None
ref_oblate = None
ns = []
reldiffsQabsprol = []
reldiffsQabspolprol = []
reldiffsQabsobl = []
reldiffsQabspolobl = []
maxdiffsQabsprol = []
maxdiffsQabspolprol = []
maxdiffsQabsobl = []
maxdiffsQabspolobl = []
for ngauss in range(nmin, nmax + 1):
    ngauss = nmax + 1 - ngauss
    output_prolate = np.fromfile(
        "ScatteringCorrection_n{0:02d}_prolate.dat".format(ngauss)
    ).reshape((len(thetas), 2))
    output_oblate = np.fromfile(
        "ScatteringCorrection_n{0:02d}_oblate.dat".format(ngauss)
    ).reshape((len(thetas), 2))

    if ref_prolate is None:
        ref_prolate = output_prolate
        ref_oblate = output_oblate
    else:
        ns.append(ngauss)
        reldiffsQabsprol.append(
            np.abs(output_prolate[:, 0] - ref_prolate[:, 0]).mean()
        )
        reldiffsQabspolprol.append(
            np.abs(output_prolate[:, 1] - ref_prolate[:, 1]).mean()
        )
        reldiffsQabsobl.append(
            np.abs(output_oblate[:, 0] - ref_oblate[:, 0]).mean()
        )
        reldiffsQabspolobl.append(
            np.abs(output_oblate[:, 1] - ref_oblate[:, 1]).mean()
        )
        maxdiffsQabsprol.append(
            np.abs(output_prolate[:, 0] - ref_prolate[:, 0]).max()
        )
        maxdiffsQabspolprol.append(
            np.abs(output_prolate[:, 1] - ref_prolate[:, 1]).max()
        )
        maxdiffsQabsobl.append(
            np.abs(output_oblate[:, 0] - ref_oblate[:, 0]).max()
        )
        maxdiffsQabspolobl.append(
            np.abs(output_oblate[:, 1] - ref_oblate[:, 1]).max()
        )

    if ngauss <= nmax_plot:
        ax[0].plot(
            thetas,
            output_prolate[:, 0],
            color=cm((ngauss - nmin) / (nmax_plot - nmin)),
        )
        ax[0].plot(
            thetas,
            output_oblate[:, 0],
            "--",
            color=cm((ngauss - nmin) / (nmax_plot - nmin)),
        )
        ax[1].plot(
            thetas,
            output_prolate[:, 1],
            color=cm((ngauss - nmin) / (nmax_plot - nmin)),
        )
        ax[1].plot(
            thetas,
            output_oblate[:, 1],
            "--",
            color=cm((ngauss - nmin) / (nmax_plot - nmin)),
        )

ax[0].plot([], [], "k-", label="prolate")
ax[0].plot([], [], "k--", label="oblate")
ax[0].legend(loc="best")

cbar = pl.colorbar(cb, cax=ax[2], label="number of quadrature points")
cbar.set_ticks(range(nmin, nmax_plot + 1))

ax[0].set_xticklabels([])
ax[1].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[1].set_xticklabels(["$0$", "$\pi{}/4$", "$\pi{}/2$"])
ax[1].set_xlabel("$\\theta{}$")
ax[0].set_ylabel("$Q_{abs}$")
ax[1].set_ylabel("$Q_{abs,pol}$")
pl.tight_layout()
pl.savefig("ScatteringCorrection.png", dpi=300, bbox_inches="tight")
pl.savefig("ScatteringCorrection.pdf", dpi=300, bbox_inches="tight")
pl.close()

pl.rcParams["figure.figsize"] = (4.5, 2.5)

pl.semilogy(
    ns,
    reldiffsQabsprol / ref_prolate[:, 0].mean(),
    "-",
    color="C0",
    label="prolate $Q_{abs}$",
)
pl.semilogy(
    ns,
    reldiffsQabspolprol / abs(ref_prolate[:, 1].mean()),
    "-",
    color="C1",
    label="prolate $Q_{abs,pol}$",
)
pl.semilogy(
    ns,
    reldiffsQabsobl / ref_oblate[:, 0].mean(),
    "-",
    color="C2",
    label="oblate $Q_{abs}$",
)
pl.semilogy(
    ns,
    reldiffsQabspolobl / abs(ref_oblate[:, 1].mean()),
    "-",
    color="C3",
    label="oblate $Q_{abs,pol}$",
)
pl.semilogy(ns, maxdiffsQabsprol / ref_prolate[:, 0].mean(), "--", color="C0")
pl.semilogy(
    ns, maxdiffsQabspolprol / abs(ref_prolate[:, 1].mean()), "--", color="C1"
)
pl.semilogy(ns, maxdiffsQabsobl / ref_oblate[:, 0].mean(), "--", color="C2")
pl.semilogy(
    ns, maxdiffsQabspolobl / abs(ref_oblate[:, 1].mean()), "--", color="C3"
)

meanline, = pl.plot([], [], "k-")
maxline, = pl.plot([], [], "k--")
leg = pl.legend([meanline, maxline], ["mean", "max"], loc="upper right")

pl.xlabel("number of quadrature points")
pl.ylabel("relative error")
pl.legend(loc="best")
pl.gca().add_artist(leg)
pl.gca().set_xticks(range(1, 20)[::2])
pl.tight_layout()
pl.savefig("ScatteringCorrection_convergence.png", dpi=300, bbox_inches="tight")
pl.savefig("ScatteringCorrection_convergence.pdf", dpi=300, bbox_inches="tight")
