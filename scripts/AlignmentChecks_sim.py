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


def run_sim(theta, shape, ad):
    return CosTuuM.get_table(
        types=CosTuuM.SILICON,
        sizes=1.0e-6,
        wavelengths=wavelengths,
        thetas=theta,
        shape_distribution=CosTuuM.SingleShapeShapeDistribution(shape),
        alignment_distribution=ad,
        dust_properties=CosTuuM.DraineDustProperties(),
        minimum_order=10,
        maximum_order=146,
        gauss_legendre_factor=2,
        tolerance=1.0e-4,
        number_of_quadrature_angles=20,
        number_of_threads=8,
        verbose=True,
        account_for_scattering=True,
        maximum_memory_size=40000000000,
    )


wavelengths = np.logspace(-5.0, -3.0, 100)

for d in [0.5, 2.0]:
    # raw
    theta = np.array([0.0, 0.5 * np.pi])
    ad = CosTuuM.SizeBasedAlignmentDistribution(
        1.0e-7, CosTuuM.DISABLE_ALIGNMENT
    )
    output = run_sim(theta, d, ad)
    output.tofile("AlignmentChecks_d{0}_raw.dat".format(d))

    # Gauss-Legendre theta quadrature angles, shared by mish and perf
    costheta, _ = np.polynomial.legendre.leggauss(50)
    theta = np.arccos(costheta)
    # mish
    ad = CosTuuM.SizeBasedAlignmentDistribution(
        1.0e-7, CosTuuM.MISHCHENKO_ALIGNMENT
    )
    output = run_sim(theta, d, ad)
    output.tofile("AlignmentChecks_d{0}_mish.dat".format(d))

    # perf
    ad = CosTuuM.SizeBasedAlignmentDistribution(
        1.0e-7, CosTuuM.DAVIS_GREENSTEIN_ALIGNMENT
    )
    output = run_sim(theta, d, ad)
    output.tofile("AlignmentChecks_d{0}_perf.dat".format(d))

    # rand
    theta = np.linspace(0.0, 0.5 * np.pi, 50)
    ad = CosTuuM.SizeBasedAlignmentDistribution(1.0, 0)
    output = run_sim(theta, d, ad)
    output.tofile("AlignmentChecks_d{0}_rand.dat".format(d))
