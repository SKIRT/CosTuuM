################################################################################
 # This file is part of CosTuuM
 # Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
import scipy.special as special


def p_exp(beta, pn):
    result = np.zeros(beta.shape)
    for i in range(len(pn)):
        result += (
            0.5 * (2.0 * i + 1.0) * pn[i, 1] * special.legendre(i)(np.cos(beta))
        )
    return result


ref = np.loadtxt("test_orientationdistribution_ref.txt")
pn = np.loadtxt("test_orientationdistribution.txt")

fig, ax = pl.subplots(2, 1, sharex=True)

xrange = ref[:, 0]
perange = p_exp(xrange, pn)
prange = ref[:, 1]
ax[0].plot(xrange, perange, "o", label="code")
ax[0].plot(xrange, prange, label="ref")
ax[0].legend(loc="best")

ax[1].semilogy(xrange, abs(perange - prange) / abs(perange + prange))

ax[1].set_xlabel("$x$")
ax[0].set_ylabel("$p(x)$")
ax[1].set_ylabel("relative difference")

pl.tight_layout()
pl.savefig("test_orientationdistribution.png", dpi=300, bbox_inches="tight")
