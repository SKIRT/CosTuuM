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
import re
import CosTuuM

nradexp = re.compile("(\d+).*NRAD")
nwavexp = re.compile("(\d+).*NWAV")
radexp = re.compile("(\d+\.\d+E[+-]\d+).*radius\(micron\)")

ad = CosTuuM.SizeBasedAlignmentDistribution(1.0e-4, 0)

# read the optical properties file
# we do several passes
# first pass: read number of radii, number of wavelengths and radii values
# we use regular expressions to parse those
# open the file for the first pass
dfile = open("1993Draine_optical_properties_C.txt")
lines = dfile.readlines()
nrad = -1
nwav = -1
rad = []
offsets = []
# loop over the lines and parse them
for i in range(len(lines)):
    line = lines[i]
    nradm = nradexp.findall(line)
    if len(nradm) > 0:
        nrad = int(nradm[0])
    nwavm = nwavexp.findall(line)
    if len(nwavm) > 0:
        nwav = int(nwavm[0])
    radm = radexp.findall(line)
    if len(radm) > 0:
        rad.append(float(radm[0]))
        offsets.append(i)
# check that we found all values we need
if nrad == -1:
    print("No number of radii found!")
    exit()
if nwav == -1:
    print("No number of wavelengths found!")
    exit()
if len(rad) != nrad:
    print("Number of radii found does not match expectation!")
    exit()

# construct the property data cube
propdata = np.zeros((nrad, nwav, 3))
# open the file again for every radius, using np.loadtxt and the
# skiprows/max_rows arguments
for i in range(len(offsets)):
    raddata = np.loadtxt(
        "1993Draine_optical_properties_C.txt",
        skiprows=offsets[i] + 2,
        max_rows=nwav,
    )
    # set the values in the data array
    propdata[i, :, 0] = rad[i] * 1.0e-6
    propdata[i, :, 1] = raddata[:, 0] * 1.0e-6
    propdata[i, :, 2] = raddata[:, 1]

dp = CosTuuM.DraineDustProperties(25.0)

isize = 10
sizes = propdata[:, 0, 0]
wavelengths = propdata[0, :, 1]
idx = wavelengths >= 1.0e-5
wavelengths = wavelengths[idx][::-1]

output_parallel = CosTuuM.get_table(
    types=CosTuuM.CARBON_PARALLEL,
    sizes=sizes[isize],
    wavelengths=wavelengths,
    thetas=0.5 * np.pi,
    shape_distribution=CosTuuM.SingleShapeShapeDistribution(1.00001),
    alignment_distribution=ad,
    dust_properties=dp,
    minimum_order=10,
    maximum_order=146,
    gauss_legendre_factor=2,
    tolerance=1.0e-4,
    number_of_threads=4,
    verbose=True,
    account_for_scattering=True,
    maximum_memory_size=5000000000,
)
output_perpendicular = CosTuuM.get_table(
    types=CosTuuM.CARBON_PERPENDICULAR,
    sizes=sizes[isize],
    wavelengths=wavelengths,
    thetas=0.5 * np.pi,
    shape_distribution=CosTuuM.SingleShapeShapeDistribution(1.00001),
    alignment_distribution=ad,
    dust_properties=dp,
    minimum_order=10,
    maximum_order=146,
    gauss_legendre_factor=2,
    tolerance=1.0e-4,
    number_of_threads=4,
    verbose=True,
    account_for_scattering=True,
    maximum_memory_size=5000000000,
)

output_parallel.tofile("DraineCarbon_parallel.dat")
output_perpendicular.tofile("DraineCarbon_perpendicular.dat")
