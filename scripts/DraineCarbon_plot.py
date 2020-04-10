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
import re

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 2.5)

nradexp = re.compile("(\d+).*NRAD")
nwavexp = re.compile("(\d+).*NWAV")
radexp = re.compile("(\d+\.\d+E[+-]\d+).*radius\(micron\)")

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

isize = 10
sizes = propdata[:, 0, 0]
wavelengths = propdata[0, :, 1]
idx = wavelengths >= 1.0e-5
wavelengths = wavelengths[idx][::-1]

output_parallel = np.fromfile("DraineCarbon_parallel.dat").reshape(
    (len(wavelengths), 2)
)
output_perpendicular = np.fromfile("DraineCarbon_perpendicular.dat").reshape(
    (len(wavelengths), 2)
)

l1, = pl.loglog(
    wavelengths * 1.0e6,
    output_parallel[:, 0],
    label="CosTuuM, $\\varepsilon{} \\parallel{} c$",
)
l2, = pl.loglog(
    wavelengths * 1.0e6,
    output_perpendicular[:, 0],
    label="CosTuuM, $\\varepsilon{} \\perp{} c$",
)
l3, = pl.loglog(
    wavelengths * 1.0e6,
    propdata[isize, idx, 2][::-1],
    label="Draine \& Lee (1984)",
)
l4, = pl.loglog(
    wavelengths * 1.0e6,
    output_parallel[:, 0] / 3.0 + 2.0 * output_perpendicular[:, 0] / 3.0,
    label="CosTuuM, $1/3-2/3$",
)
pl.xlim(10.0, 1000.0)
pl.xlabel("$\\lambda{}$ ($\\mu{}$m)")
pl.ylabel("$Q_{abs}$")
first_legend = pl.legend(handles=[l1, l2], loc="upper right")
pl.gca().add_artist(first_legend)
pl.legend(handles=[l3, l4], loc="lower left")
pl.title("$a = {0:.1f}~\\mu{{}}$m".format(sizes[isize] * 1.0e6))
pl.tight_layout()
pl.savefig("1993Draine_carbon.png", dpi=300, bbox_inches="tight")
pl.savefig("1993Draine_carbon.pdf", dpi=300, bbox_inches="tight")
