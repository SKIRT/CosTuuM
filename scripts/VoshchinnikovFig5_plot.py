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

vosh = np.load("1993Voshchinnikov_fig5.npz", mmap_mode="r")

mr = 1.31 + 0.01j
ds = vosh["d"]
thetas = vosh["theta"]
Qvosh = vosh["Qs"]
rV = 1.5e-7
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

cost = np.load("VoshchinnikovFig5.npz", mmap_mode="r")
Qcos = cost["Qs"]

colors = ["C0", "C1"]

fig, ax = pl.subplots(2, 2, sharex="col", sharey="row")
ax[0][0].set_title("prolate")
for itheta in range(len(thetas)):
    ax[0][0].plot(
        1.0 / ds[:50],
        Qcos[:50, itheta, 0],
        "-",
        color=colors[itheta],
        alpha=0.5,
    )
    ax[0][0].plot(
        1.0 / ds[:50:4], Qvosh[:50:4, itheta, 0], "x", color=colors[itheta]
    )
    ax[1][0].plot(
        1.0 / ds[:50],
        Qcos[:50, itheta, 1],
        "-",
        color=colors[itheta],
        alpha=0.5,
    )
    ax[1][0].plot(
        1.0 / ds[:50:4], Qvosh[:50:4, itheta, 1], "x", color=colors[itheta]
    )
ax[0][1].set_title("oblate")
for itheta in range(len(thetas)):
    ax[0][1].plot(
        ds[50:], Qcos[50:, itheta, 0], "-", color=colors[itheta], alpha=0.5
    )
    ax[0][1].plot(ds[50::4], Qvosh[50::4, itheta, 0], "x", color=colors[itheta])
    ax[1][1].plot(
        ds[50:],
        np.abs(Qcos[50:, itheta, 1]),
        "-",
        color=colors[itheta],
        alpha=0.5,
    )
    ax[1][1].plot(
        ds[50::4], np.abs(Qvosh[50::4, itheta, 1]), "x", color=colors[itheta]
    )

ax[0][0].plot([], [], "-", color=colors[0], label="$\\theta{}=0$")
ax[0][0].plot([], [], "-", color=colors[1], label="$\\theta{}=\\pi{}/2$")
ax[0][0].legend(loc="best", handlelength=0.5, labelspacing=0.1)

ax[1][0].plot([], [], "k-", label="\\textsc{CosTuuM}")
ax[1][0].legend(loc="best", handlelength=0.5, labelspacing=0.1)
ax[1][1].plot([], [], "kx", label="V \& F (1993)")
ax[1][1].legend(loc="best", handlelength=0.5, labelspacing=0.1)

ax[0][0].set_ylabel("$Q_{ext}$")
ax[1][0].set_ylabel("$|Q_{ext,pol}|$")
ax[1][0].set_xlabel("$1/d$")
ax[1][1].set_xlabel("$d$")

pl.tight_layout()
pl.savefig("VoshchinnikovFig5.pdf", dpi=300, bbox_inches="tight")
