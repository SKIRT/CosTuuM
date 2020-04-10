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


def get_L(d):
    if d == 0.0:
        return 0.0
    elif d < 1.0:
        e = np.sqrt(1.0 - d ** 2)
        return (
            (1.0 - e ** 2)
            / e ** 2
            * (0.5 * np.log((1.0 + e) / (1.0 - e)) / e - 1.0)
        )
    elif d > 1.0:
        e = np.sqrt(1.0 - 1.0 / d ** 2)
        return (1.0 - np.sqrt(1.0 - e ** 2) / e * np.arcsin(e)) / e ** 2
    else:
        return 1.0 / 3.0


def dLdd(drange):
    erange = np.sqrt(
        np.where(drange < 1.0, 1.0 - drange ** 2, 1.0 - 1.0 / drange ** 2)
    )
    Prange = np.where(
        drange < 1.0,
        0.5
        * np.sqrt(1.0 - erange ** 2)
        / erange ** 5
        * (
            (3.0 - erange ** 2) * np.log((1.0 + erange) / (1.0 - erange))
            - 6.0 * erange
        ),
        (1.0 - erange ** 2)
        / erange ** 5
        * (
            (3.0 - 2.0 * erange ** 2) * np.arcsin(erange)
            - 3.0 * erange * np.sqrt(1.0 - erange ** 2)
        ),
    )
    return np.where(drange == 0.0, np.zeros(drange.shape), Prange)


def cde2(L):
    return 12.0 * L * (1.0 - L) ** 2


drange = np.linspace(0.0, 7.5, 1000)
ifilt = (drange > 0.19) & (drange < 6.96)
Lrange = np.array([get_L(d) for d in drange])
dLddP = dLdd(drange)
cde2P = cde2(Lrange) * dLddP
print((drange * cde2P).sum() * ((drange[-1] - drange[0]) / len(drange)))
print(cde2P[drange < 1.0].sum() * (drange[1] - drange[0]))
print(dLddP[drange < 1.0].sum() * (drange[1] - drange[0]))
pl.plot(drange, cde2P, label="CDE2")
pl.fill_between(
    drange[ifilt], np.zeros(drange[ifilt].shape), cde2P[ifilt], alpha=0.5
)
pl.plot(drange, dLddP, label="CDS")
pl.gca().axvline(x=1.0, linestyle="--", color="k")
pl.text(-0.3, 0.55, "prolate")
pl.text(1.2, 0.5, "oblate")
pl.legend(loc="best")
pl.xlabel("$d = a/b$")
pl.ylabel("$P(d)$")
pl.tight_layout()
pl.savefig("shapes.png", dpi=300, bbox_inches="tight")
pl.savefig("shapes.pdf", dpi=300, bbox_inches="tight")
