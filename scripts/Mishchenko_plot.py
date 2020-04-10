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
pl.rcParams["figure.figsize"] = (9, 5)

nice_names = {
    "pdg": "Perfect Davis-Greenstein alignment",
    "ia": "Imperfect Mishchenko alignment",
    "ro": "Random orientations (no alignment)",
}
axis_ratios = {"prolate": 0.5, "oblate": 2.0}
orientations = {
    "pdg": {"prolate": 0.0, "oblate": 1.0},
    "ia": {"prolate": 1.0 / 5.0, "oblate": 3.0 / 5.0},
    "ro": {"prolate": 1.0 / 3.0, "oblate": 1.0 / 3.0},
}
r_ev = 2.0e-7
thetas = np.array([0.0, np.pi / 6.0, np.pi / 3.0, np.pi / 2.0])
finethetas = np.linspace(0.0, np.pi / 2.0, 100)
wavelengths = np.array([0.2, 0.55, 1.0, 5.0]) * 1.0e-6
colors = ["C0", "C1", "C2", "C3"]


def plot_case(axis_ratio_key, orientation_key, ax0, ax1, symbol, line):
    fname = "Mishchenko_{0}_{1}".format(axis_ratio_key, orientation_key)
    refdata = np.loadtxt("1991" + fname + ".txt")
    mdata = np.zeros((4, 4, 3))
    mdata[:, 0, :] = refdata[:4, 1:]
    if orientation_key == "ro":
        mdata[:, 1, :] = refdata[:4, 1:]
        mdata[:, 2, :] = refdata[:4, 1:]
        mdata[:, 3, :] = refdata[:4, 1:]
    else:
        mdata[:, 1, :] = refdata[4:8, 1:]
        mdata[:, 2, :] = refdata[8:12, 1:]
        mdata[:, 3, :] = refdata[12:, 1:]

    results = np.fromfile(fname + ".dat").reshape(
        (len(wavelengths), len(finethetas), 3)
    )

    for i in range(len(wavelengths)):
        ax0.plot(thetas, mdata[i, :, 0], "o", color=colors[i], marker=symbol)
        ax0.plot(finethetas, results[i, :, 0], color=colors[i], linestyle=line)
        ax1.plot(thetas, mdata[i, :, 1], "o", color=colors[i], marker=symbol)
        ax1.plot(finethetas, results[i, :, 1], color=colors[i], linestyle=line)


fig, ax = pl.subplots(2, 3, sharey="row", sharex=True)
symbols = ["x", "."]
lines = ["-", ":"]
for iar in range(len(axis_ratios.keys())):
    axis_ratio_key = [*axis_ratios][iar]
    for ior in range(len(orientations.keys())):
        orientation_key = [*orientations][ior]
        plot_case(
            axis_ratio_key,
            orientation_key,
            ax[0][ior],
            ax[1][ior],
            symbols[iar],
            lines[iar],
        )

for ior in range(len(orientations.keys())):
    orientation_key = [*orientations][ior]
    ax[0][ior].set_title(nice_names[orientation_key])
ax[0][0].set_ylabel("$Q_{ext}$")
ax[1][0].set_ylabel("$Q_{ext,pol}$")
ax[1][0].set_xlabel("$\\theta{}$")
ax[1][1].set_xlabel("$\\theta{}$")
ax[1][2].set_xlabel("$\\theta{}$")
ax[1][0].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[1][0].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])

for i in range(len(wavelengths)):
    ax[1][2].plot(
        [],
        [],
        "o",
        color=colors[i],
        label="$\\lambda{{}}={0:.2f}~\\mu{{}}$m".format(wavelengths[i] * 1.0e6),
    )
for iar in range(len(axis_ratios.keys())):
    axis_ratio_key = [*axis_ratios][iar]
    ax[1][2].plot(
        [],
        [],
        color="k",
        marker=symbols[iar],
        linestyle=lines[iar],
        label="$d={0:.1f}$".format(axis_ratios[axis_ratio_key]),
    )
ax[1][2].legend(loc="best")

pl.tight_layout()
pl.savefig("1991Mishchenko.png", dpi=300, bbox_inches="tight")
pl.savefig("1991Mishchenko.pdf", dpi=300, bbox_inches="tight")
