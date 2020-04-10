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
import scipy.interpolate as interpol
import scipy.integrate as integ
import re

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 2.5)

# integrand for the integral in the expression for the active surface area
# ratio
def active_surface_area_integrand(x, d4m1, d2m1, d):
    return np.sqrt((d4m1 * x ** 2 + 1.0) / (d2m1 * x ** 2 + 1.0)) / np.cbrt(d)


# compute the ratio of the actual average active surface area of a spheroid
# with axis ratio d, and the naive average active surface area pi*rV**2
def get_active_surface_area_factor(d):
    d4m1 = d ** 4 - 1.0
    d2m1 = d ** 2 - 1.0
    return (
        0.5
        * integ.quad(
            active_surface_area_integrand, -1.0, 1.0, args=(d4m1, d2m1, d)
        )[0]
    )


# read the 1993Draine_optical_properties.txt file
# we do several passes
# first pass: read number of radii, number of wavelengths and radii values
# we use regular expressions to parse those
nradexp = re.compile("(\d+).*NRAD")
nwavexp = re.compile("(\d+).*NWAV")
radexp = re.compile("(\d+\.\d+E[+-]\d+).*radius\(micron\)")
# open the file for the first pass
dfile = open("1993Draine_optical_properties.txt")
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
        "1993Draine_optical_properties.txt",
        skiprows=offsets[i] + 2,
        max_rows=nwav,
    )
    # set the values in the data array
    propdata[i, :, 0] = rad[i] * 1.0e-6
    propdata[i, :, 1] = raddata[:, 0] * 1.0e-6
    propdata[i, :, 2] = raddata[:, 1]

# parse the dielectric function file
epsdata = np.loadtxt("1993Draine_dielectric_function.txt", skiprows=6)
epsdata[:, 3] += 1.0

# construct a linear interpolation function
mr = interpol.interp1d(
    epsdata[:, 0] * 1.0e-6, epsdata[:, 3] + epsdata[:, 4] * 1.0j
)

# create a constant refractive index function
def get_refractive_index(wavelength, grain_size, grain_type, interpfunc):
    return interpfunc(5.0e-4)


wavelengths = np.unique(propdata[0, :, 1])
wavelengths = wavelengths[wavelengths >= 1.0e-5]

size_index = 15
sizes = np.unique(propdata[size_index, 0, 0])

d1 = np.fromfile("Draine_mrconst_d1.dat").reshape((len(wavelengths), 2))
d05 = np.fromfile("Draine_mrconst_d05.dat").reshape((len(wavelengths), 2))
d2 = np.fromfile("Draine_mrconst_d2.dat").reshape((len(wavelengths), 2))

area_correction = 1.0 / get_active_surface_area_factor(0.5)
pl.semilogx(
    wavelengths * 1.0e6,
    d05[:, 0] / d1[:, 0] * area_correction,
    "--",
    color="C0",
    label="$d=1/2$",
)
area_correction = 1.0 / get_active_surface_area_factor(2.0)
pl.semilogx(
    wavelengths * 1.0e6,
    d2[:, 0] / d1[:, 0] * area_correction,
    "--",
    color="C1",
    label="$d=2$",
)
pl.gca().axhline(y=1.0, linestyle="--", color="k")
pl.xlim(10.0, 1000.0)
pl.xlabel("$\\lambda{}$ ($\\mu{}$m)")
pl.ylabel("$C_A Q_{abs}/Q_{abs,sph}$")
pl.legend(loc="best")
refractive_index = get_refractive_index(0.0, 0.0, 0.0, mr)
pl.title(
    "$a = {0:.1f}~\\mu{{}}$m, $m_r={1:.2f}+{2:.2f}i$".format(
        sizes[0] * 1.0e6, refractive_index.real, refractive_index.imag
    )
)
pl.tight_layout()
pl.savefig("1993Draine_correction.png", dpi=300, bbox_inches="tight")
pl.savefig("1993Draine_correction.pdf", dpi=300, bbox_inches="tight")
