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
import CosTuuM

vosh = np.load("1993Voshchinnikov_fig5.npz", mmap_mode="r")

mr = 1.31 + 0.01j
ds = vosh["d"]
thetas = vosh["theta"]
Qs = vosh["Qs"]
rV = 1.5e-7
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

Qs = np.zeros((ds.shape[0], thetas.shape[0], 2))
for id in range(len(ds)):
    d = ds[id]
    print("d:", d)

    def get_refractive_index(wavelength, size, gtype):
        return np.array(mr)

    dust_properties = CosTuuM.CustomDustProperties(get_refractive_index)
    results = CosTuuM.get_table(
        CosTuuM.SILICON,
        rV,
        wav,
        thetas,
        CosTuuM.SingleShapeShapeDistribution(d),
        CosTuuM.SizeBasedAlignmentDistribution(
            rV - 10.0, CosTuuM.DISABLE_ALIGNMENT
        ),
        dust_properties,
        do_extinction=True,
        do_absorption=False,
    )

    Qs[id, :, 0] = results[:, 0]
    Qs[id, :, 1] = results[:, 1]


np.savez("VoshchinnikovFig5.npz", theta=thetas, Qs=Qs)
