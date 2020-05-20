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

vosh = np.load("1993Voshchinnikov_fig4.npz", mmap_mode="r")

mrs = vosh["mr"]
ds = vosh["d"]
thetas = vosh["theta"]
Qvosh = vosh["Qpol"]
rV = 1.5e-7
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

Qcos = np.zeros((mrs.shape[0], ds.shape[0], thetas.shape[0]))
for imr in range(len(mrs)):
    mr = mrs[imr]
    print("mr:", mr)
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

        Qcos[imr, id, :] = results[:, 1]


np.savez("VoshchinnikovFig4.npz", theta=thetas, Qcos=Qcos)
