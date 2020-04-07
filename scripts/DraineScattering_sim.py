import numpy as np
import scipy.interpolate as interpol
import re
import CosTuuM

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

# construct a custom dust properties object that uses the interpolation function
def get_refractive_index(wavelength, grain_size, grain_type, interpfunc):
    return interpfunc(wavelength)


dp = CosTuuM.CustomDustProperties(get_refractive_index, mr)

ad = CosTuuM.SizeBasedAlignmentDistribution(1.0e-4, 0)

size_index = 20
sizes = np.unique(propdata[size_index, 0, 0])
wavelengths = np.unique(propdata[0, :, 1])
wavelengths = wavelengths[wavelengths >= 1.0e-5]

without_scattering = CosTuuM.get_table(
    types=CosTuuM.SILICON,
    sizes=sizes,
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
    account_for_scattering=False,
    maximum_memory_size=5000000000,
)
with_scattering = CosTuuM.get_table(
    types=CosTuuM.SILICON,
    sizes=sizes,
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

without_scattering.tofile("Draine_without_scattering.dat")
with_scattering.tofile("Draine_with_scattering.dat")