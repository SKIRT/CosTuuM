import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import scipy.interpolate as interpol
import scipy.integrate as integ
import re
import CosTuuM

pl.rcParams["text.usetex"] = True

# integrand for the integral in the expression for the active surface area
# ratio
def active_surface_area_integrand(x, d4m1, d2m1):
    return np.sqrt((d4m1 * x ** 2 + 1.0) / (d2m1 * x ** 2 + 1.0))


# compute the ratio of the actual average active surface area of a spheroid
# with axis ratio d, and the naive average active surface area pi*rV**2
def get_active_surface_area_factor(d):
    d4m1 = d ** 4 - 1.0
    d2m1 = d ** 2 - 1.0
    return (
        0.5
        * integ.quad(
            active_surface_area_integrand, -1.0, 1.0, args=(d4m1, d2m1)
        )[0]
        / np.cbrt(d)
    )


# get the active surface area of the spheroid with the given equal volume
# sphere radius and axis ratio
def get_active_surface_area(rV, d):
    return np.pi * rV ** 2 * get_active_surface_area_factor(d)


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
# construct a linear interpolation function
mr = interpol.interp1d(
    epsdata[:, 0] * 1.0e-6, 1.0 + epsdata[:, 3] + epsdata[:, 4] * 1.0j
)

lcol = 10
wcol = 40

ai = propdata[lcol, wcol, 0]
wi = propdata[lcol, wcol, 1]
mri = mr(wi)

thetas = np.linspace(0.0, np.pi, 20)

Tmatrix = CosTuuM.TMatrix(
    particle_radius=ai,
    axis_ratio=1.00001,
    wavelength=wi,
    refractive_index=mri,
    cos2beta=1.0,
)
Qabsd1 = Tmatrix.get_average_absorption_cross_section() / (np.pi * ai ** 2)
Qabsd1theta = Tmatrix.get_absorption_cross_sections(theta=thetas) / (
    np.pi * ai ** 2
)
Tmatrix = CosTuuM.TMatrix(
    particle_radius=ai,
    axis_ratio=2.0,
    wavelength=wi,
    refractive_index=mri,
    cos2beta=1.0,
    is_equal_volume_radius=False,
)
Qabsd2 = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
    ai, 2.0
)
Qabsd2theta = Tmatrix.get_absorption_cross_sections(
    theta=thetas
) / get_active_surface_area(ai, 2.0)
Tmatrix = CosTuuM.TMatrix(
    particle_radius=ai,
    axis_ratio=2.0,
    wavelength=wi,
    refractive_index=mri,
    cos2beta=3.0 / 5.0,
    is_equal_volume_radius=False,
)
Qabsd2ia = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
    ai, 2.0
)
Qabsd2thetaia = Tmatrix.get_absorption_cross_sections(
    theta=thetas
) / get_active_surface_area(ai, 2.0)
Tmatrix = CosTuuM.TMatrix(
    particle_radius=ai,
    axis_ratio=0.5,
    wavelength=wi,
    refractive_index=mri,
    cos2beta=1.0,
    is_equal_volume_radius=False,
)
Qabsd05 = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
    ai, 0.5
)
Qabsd05theta = Tmatrix.get_absorption_cross_sections(
    theta=thetas
) / get_active_surface_area(ai, 0.5)
Tmatrix = CosTuuM.TMatrix(
    particle_radius=ai,
    axis_ratio=0.5,
    wavelength=wi,
    refractive_index=mri,
    cos2beta=1.0 / 5.0,
    is_equal_volume_radius=False,
)
Qabsd05ia = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
    ai, 0.5
)
Qabsd05thetaia = Tmatrix.get_absorption_cross_sections(
    theta=thetas
) / get_active_surface_area(ai, 0.5)

fig, ax = pl.subplots(2, 1, sharex=True)

ax[0].axhline(y=propdata[lcol, wcol, 2], label="Laor \& Draine (1993)")
ax[0].axhline(y=Qabsd2, color="C2", linestyle="--")
ax[0].axhline(y=Qabsd05, color="C3", linestyle="--")
ax[0].axhline(y=Qabsd2ia, color="C4", linestyle="--")
ax[0].axhline(y=Qabsd05ia, color="C5", linestyle="--")
ax[0].plot(
    thetas,
    Qabsd1theta[:, 0],
    "o",
    color="C1",
    label="CosTuuM spherical (no alignment)",
)
ax[0].plot(
    thetas,
    Qabsd2theta[:, 0],
    "o",
    color="C2",
    label="CosTuuM oblate, perfect alignment",
)
ax[0].plot(
    thetas,
    Qabsd05theta[:, 0],
    "o",
    color="C3",
    label="CosTuuM prolate, perfect alignment",
)
ax[0].plot(
    thetas,
    Qabsd2thetaia[:, 0],
    "o",
    color="C4",
    label="CosTuuM oblate, imperfect alignment",
)
ax[0].plot(
    thetas,
    Qabsd05thetaia[:, 0],
    "o",
    color="C5",
    label="CosTuuM prolate, imperfect alignment",
)

ax[1].plot(thetas, Qabsd1theta[:, 1], "o", color="C1")
ax[1].plot(thetas, Qabsd2theta[:, 1], "o", color="C2")
ax[1].plot(thetas, Qabsd05theta[:, 1], "o", color="C3")
ax[1].plot(thetas, Qabsd2thetaia[:, 1], "o", color="C4")
ax[1].plot(thetas, Qabsd05thetaia[:, 1], "o", color="C5")

ax[0].set_title(
    "grain size: ${0:.2f}\\mu{{}}$m, $\\lambda{{}} = {1:.2f}\\mu{{}}$m".format(
        ai * 1.0e6, wi * 1.0e6
    )
)
ax[0].set_ylabel("$Q_{abs}$")
ax[1].set_xlabel("$\\theta{}$")
ax[1].set_ylabel("$Q_{abspol}$")
ax[0].legend(loc="upper right")
pl.tight_layout()
pl.savefig("1993Draine_theta.png", dpi=300)
pl.close()

lcol = 18
maxl = 100
Qabsd1 = np.zeros(maxl)
Qabsd2 = np.zeros(maxl)
Qabsd05 = np.zeros(maxl)
for iw in range(maxl):
    print(iw)
    ai = propdata[lcol, iw, 0]
    wi = propdata[lcol, iw, 1]
    mri = mr(wi)

    # using an axis ratio of exactly 1 causes small (~1%) deviations for
    # some wavelengths. We follow Mishchenko's suggestion and use a value
    # slightly higher than 1
    Tmatrix = CosTuuM.TMatrix(
        particle_radius=ai,
        axis_ratio=1.00001,
        wavelength=wi,
        refractive_index=mri,
        cos2beta=1.0 / 3.0,
    )
    Qabsd1[iw] = Tmatrix.get_average_absorption_cross_section() / (
        np.pi * ai ** 2
    )
    Tmatrix = CosTuuM.TMatrix(
        particle_radius=ai,
        axis_ratio=2.0,
        wavelength=wi,
        refractive_index=mri,
        cos2beta=1.0 / 3.0,
        is_equal_volume_radius=False,
    )
    Qabsd2[
        iw
    ] = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
        ai, 0.5
    )
    Tmatrix = CosTuuM.TMatrix(
        particle_radius=ai,
        axis_ratio=0.5,
        wavelength=wi,
        refractive_index=mri,
        cos2beta=1.0 / 3.0,
        is_equal_volume_radius=False,
    )
    Qabsd05[
        iw
    ] = Tmatrix.get_average_absorption_cross_section() / get_active_surface_area(
        ai, 0.5
    )

fig, ax = pl.subplots(2, 1, sharex=True)

ax[0].loglog(
    propdata[0, :maxl, 1] * 1.0e6, Qabsd1, "-", label="CosTuuM spherical"
)
ax[0].loglog(propdata[0, :maxl, 1] * 1.0e6, Qabsd2, "-", label="CosTuuM oblate")
ax[0].loglog(
    propdata[0, :maxl, 1] * 1.0e6, Qabsd05, "-", label="CosTuuM prolate"
)
ax[0].loglog(
    propdata[0, :maxl, 1] * 1.0e6,
    propdata[lcol, :maxl, 2],
    "-",
    label="Laor \& Draine (1993)",
)


def reldiff(a, b):
    a = np.array(a)
    b = np.array(b)
    return 2.0 * np.abs((a - b) / (a + b))


ax[1].loglog(
    propdata[lcol, :maxl, 1] * 1.0e6, reldiff(Qabsd1, propdata[lcol, :maxl, 2])
)
ax[1].loglog(
    propdata[lcol, :maxl, 1] * 1.0e6, reldiff(Qabsd2, propdata[lcol, :maxl, 2])
)
ax[1].loglog(
    propdata[lcol, :maxl, 1] * 1.0e6, reldiff(Qabsd05, propdata[lcol, :maxl, 2])
)

ax[0].set_title("grain size: ${0:.2f}\\mu{{}}$m".format(rad[lcol]))
ax[0].set_ylabel("$Q_{abs}$")
ax[1].set_xlabel("$\\lambda{}$ ($\\mu{}$m)")
ax[1].set_ylabel("reldiff")

ax[0].legend(loc="best")
pl.tight_layout()
pl.savefig("1993Draine.png", dpi=300)
