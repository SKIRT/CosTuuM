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

ad = CosTuuM.SizeBasedAlignmentDistribution(1.0e-7, 0)
dp = CosTuuM.DraineDustProperties()
thetas = np.linspace(0.0, 0.5 * np.pi, 50)

sd = CosTuuM.DraineHensleyShapeDistribution(20, 0.0, 0.96)
lims = sd.get_limits()
dprol = lims[0]
dobl = lims[1]
print(lims)

nmin = 1
nmax = 20

for ngauss in range(nmin, nmax + 1):
    output_prolate = CosTuuM.get_table(
        types=CosTuuM.SILICON,
        sizes=1.0e-6,
        wavelengths=1.0e-5,
        thetas=thetas,
        shape_distribution=CosTuuM.SingleShapeShapeDistribution(dprol),
        alignment_distribution=ad,
        dust_properties=dp,
        minimum_order=10,
        maximum_order=146,
        gauss_legendre_factor=2,
        tolerance=1.0e-4,
        number_of_quadrature_angles=ngauss,
        number_of_threads=1,
        verbose=True,
        account_for_scattering=True,
        maximum_memory_size=5000000000,
    )
    output_oblate = CosTuuM.get_table(
        types=CosTuuM.SILICON,
        sizes=1.0e-6,
        wavelengths=1.0e-5,
        thetas=thetas,
        shape_distribution=CosTuuM.SingleShapeShapeDistribution(dobl),
        alignment_distribution=ad,
        dust_properties=dp,
        minimum_order=10,
        maximum_order=146,
        gauss_legendre_factor=2,
        tolerance=1.0e-4,
        number_of_quadrature_angles=ngauss,
        number_of_threads=1,
        verbose=True,
        account_for_scattering=True,
        maximum_memory_size=5000000000,
    )

    output_prolate.tofile(
        "ScatteringCorrection_n{0:02d}_prolate.dat".format(ngauss)
    )
    output_oblate.tofile(
        "ScatteringCorrection_n{0:02d}_oblate.dat".format(ngauss)
    )
