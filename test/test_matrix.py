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

##
# @file test_matrix.py
#
# @brief Script to read and analyse the output from the testMatrix unit test.
#
# We first reconstruct the output array from the binary dump file and then
# compare it with the same matrix inverted using numpy.linalg.inv.
#
# @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
##

# load numpy for file input and matrix inversion
import numpy as np

# read the metadata from the dump file (matrix and element sizes)
metadata = np.fromfile("test_matrix.dump", dtype=np.uint32, count=3)
num_x = metadata[0]
num_y = metadata[1]

# read the actual matrix and reshape it to the correct shape
data = np.fromfile("test_matrix.dump", dtype=np.cdouble, offset=12)
data = data.reshape((num_x, num_y))

# construct the same initial matrix for the reference solution
A = np.array(
    [
        [1.0 + 1.0j, 2.0 + 1.0j, 3.0 + 1.0j],
        [4.0 + 2.0j, 5.0 + 1.0j, 6.0 + 1.0j],
        [8.0 + 3.0j, 9.0 + 1.0j, -2.0j],
    ],
    dtype=np.cdouble,
)
# invert it using numpy.linalg.inv
Ainv = np.linalg.inv(A)

# compare the unit test result and the numpy library result
for i in range(3):
    for j in range(3):
        reldiff = abs(data[i, j] - Ainv[i, j]) / abs(data[i, j] + Ainv[i, j])
        print("data:", data[i, j], "ref:", Ainv[i, j], "reldiff:", reldiff)
