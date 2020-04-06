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

refs = {
    "Astronomical Silicate": {
        "Qabsfile": "1993Draine_optical_properties.txt",
        "epsfile": "1993Draine_dielectric_function.txt",
        "rowskip": 6,
    },
    "Smoothed UV Astronomical Silicate": {
        "Qabsfile": "1993Draine_optical_properties_suvSil.txt",
        "epsfile": "1993Draine_dielectric_function_suvSil.txt",
        "rowskip": 9,
    },
    "Silicon Carbide": {
        "Qabsfile": "1993Draine_optical_properties_SiC.txt",
        "epsfile": "1993Draine_dielectric_function_SiC.txt",
        "rowskip": 6,
    },
}

for ref in refs:
    # read the optical properties file
    # we do several passes
    # first pass: read number of radii, number of wavelengths and radii values
    # we use regular expressions to parse those
    # open the file for the first pass
    dfile = open(refs[ref]["Qabsfile"])
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
            refs[ref]["Qabsfile"], skiprows=offsets[i] + 2, max_rows=nwav
        )
        # set the values in the data array
        propdata[i, :, 0] = rad[i] * 1.0e-6
        propdata[i, :, 1] = raddata[:, 0] * 1.0e-6
        propdata[i, :, 2] = raddata[:, 1]

    sizes = propdata[:, 0, 0]
    wavelengths = propdata[0, :, 1]
    idx = wavelengths >= 1.0e-5
    wavelengths = wavelengths[idx][::-1]

    reldiffmax = np.fromfile("DraineAll_" + ref.replace(" ", "_") + ".dat")

    refs[ref]["sizes"] = sizes
    refs[ref]["reldiffmax"] = reldiffmax

for ref in refs:
    pl.loglog(
        refs[ref]["sizes"] * 1.0e6, refs[ref]["reldiffmax"], "-", label=ref
    )

pl.xlabel("$a$ ($\\mu{}$m)")
pl.ylabel("maximum relative error")
pl.legend(loc="best")
pl.savefig("1993Draine_all.png", dpi=300, bbox_inches="tight")
pl.savefig("1993Draine_all.pdf", dpi=300, bbox_inches="tight")
