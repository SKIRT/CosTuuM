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
pl.rcParams["figure.figsize"] = (9.0, 4.5)

fig, ax = pl.subplots(2, 2, sharex=True, sharey=True)

wavelengths = np.logspace(-5.0, -3.0, 100) * 1.0e6
costheta, _ = np.polynomial.legendre.leggauss(50)
theta = np.arccos(costheta)
ds = [0.5, 2.0]
lambdas = [0, 65]
cos2beta = np.array([0.2, 0.6])
Rs = 2.0 * (1.5 * cos2beta - 0.5)
for id in range(len(ds)):
    for ilambda in range(len(lambdas)):
        iilambda = lambdas[ilambda]
        d = ds[id]
        R = Rs[id]
        Qraw = np.fromfile("AlignmentChecks_d{0}_raw.dat".format(d)).reshape(
            (len(wavelengths), 2, 2)
        )
        Qpar = Qraw[iilambda, 0, 0]
        Qperp = Qraw[iilambda, 1, 0]

        Qabs = (2.0 * Qperp + Qpar) / 3.0 + R * (Qpar - Qperp) * (
            1.0 / 3.0 - 0.5 * np.sin(theta) ** 2
        )
        Qabspol = -0.5 * R * (Qpar - Qperp) * np.sin(theta) ** 2

        Qmish = np.fromfile("AlignmentChecks_d{0}_mish.dat".format(d)).reshape(
            (len(wavelengths), len(theta), 2)
        )
        Qmishquad = Qmish[iilambda]

        Qrand = np.fromfile("AlignmentChecks_d{0}_rand.dat".format(d)).reshape(
            (len(wavelengths), len(theta), 2)
        )
        Qperf = np.fromfile("AlignmentChecks_d{0}_perf.dat".format(d)).reshape(
            (len(wavelengths), len(theta), 2)
        )
        Qperfquad = Qperf[iilambda]

        label = None
        if ilambda == 0 and id == 0:
            label = "Picket fence"
        ax[ilambda][id].plot(theta, np.abs(Qabspol) / Qabs, "-", label=label)
        for f in [0.2, 0.4, 0.6]:
            Qmix = (1.0 - f) * Qrand[iilambda, 0, :] + f * Qperfquad
            label = None
            if ilambda == 0 and id == 1:
                label = "$f_A={0}$".format(f)
            ax[ilambda][id].plot(
                theta, np.abs(Qmix[:, 1]) / Qmix[:, 0], "-", label=label
            )
        label = None
        if ilambda == 0 and id == 0:
            label = "CosTuuM"
        t = 1
        ax[ilambda][id].plot(
            theta[::t],
            np.abs(Qmishquad[::t, 1]) / Qmishquad[::t, 0],
            ".",
            label=label,
        )

    ax[0][id].set_title("$d={0}$".format(d))

ax[0][0].set_ylabel("$\\langle{}Q_{\\rm{}abs}\\rangle{}$")
for ilambda in range(len(lambdas)):
    ax[ilambda][0].set_ylabel(
        "$P_L(\\lambda{{}}={0:.0f}~\\mu{{}}$m$)$".format(
            wavelengths[lambdas[ilambda]]
        )
    )
ax[1][0].set_xlabel("$\\theta{}$")
ax[1][1].set_xlabel("$\\theta{}$")
ax[0][0].legend(loc="upper left")
ax[0][1].legend(
    loc="upper left", title="Dyck \& Beichman", ncol=3, columnspacing=1.0
)
ax[1][0].set_xlim(0.0, 0.5 * np.pi)
# ax[1][0].set_ylim(0.0, 0.38)
ax[1][0].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[1][0].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])
pl.tight_layout()
pl.savefig("AlignmentChecks_zenith.png", dpi=300, bbox_inches="tight")
pl.savefig("AlignmentChecks_zenith.pdf", dpi=300, bbox_inches="tight")
