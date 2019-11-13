import numpy as np
import CTMmodule
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True

mrdata = np.loadtxt("1985Draine.txt")
mrdict = {}
for row in mrdata:
    nnorm = np.sqrt(row[1] ** 2 + row[2] ** 2)
    n = np.sqrt(0.5 * (nnorm + row[1]))
    k = np.sqrt(0.5 * (nnorm - row[1]))
    mrdict[row[0]] = n + 1.0j * k

fname = "1991Mishchenko_oblate_ia.txt"

rawdata = np.loadtxt(fname)

data = np.zeros((rawdata.shape[0], rawdata.shape[1] + 1))
data[:, 1:6] = rawdata[:, :]
data[:4, 0] = 0.0
data[4:8, 0] = np.pi / 6.0
data[8:12, 0] = np.pi / 3.0
data[12:16, 0] = 0.5 * np.pi

rdict = {}
for row in mrdata:
    rdict[row[0]] = {
        "theta": [],
        "Qext": [],
        "Qextref": [],
        "Qpol": [],
        "Qpolref": [],
        "Qcpol": [],
        "Qcpolref": [],
    }
for row in data:
    Tmatrix = CTMmodule.TMatrix(
        particle_radius=2.0e-7,
        axis_ratio=2.0,
        wavelength=row[1] * 1.0e-6,
        refractive_index=mrdict[row[1]],
    )
    K = Tmatrix.get_extinction_matrix(theta=row[0], phi=0.0)

    parea = np.pi * 2.0e-7 ** 2
    Qext = K[0, 0] / parea
    Qpol = K[0, 1] / parea
    Qcpol = K[2, 3] / parea

    rdict[row[1]]["theta"].append(row[0])
    rdict[row[1]]["Qext"].append(Qext)
    rdict[row[1]]["Qextref"].append(row[2])
    rdict[row[1]]["Qpol"].append(Qpol)
    rdict[row[1]]["Qpolref"].append(row[3] / 10.0)
    rdict[row[1]]["Qcpol"].append(Qcpol)
    rdict[row[1]]["Qcpolref"].append(row[4] / 10.0)

fig, ax = pl.subplots(3, 2, sharex=True)


def reldiff(a, b):
    a = np.array(a)
    b = np.array(b)
    return 2.0 * np.abs((a - b) / (a + b))


ccount = 0
for l in rdict:
    ax[0][0].plot(
        rdict[l]["theta"], rdict[l]["Qext"], "o", color="C%i" % ccount
    )
    ax[0][0].plot(
        rdict[l]["theta"], rdict[l]["Qextref"], "-", color="C%i" % ccount
    )
    ax[0][1].semilogy(
        rdict[l]["theta"],
        reldiff(rdict[l]["Qext"], rdict[l]["Qextref"]),
        "-",
        color="C%i" % ccount,
    )
    ax[1][0].plot(
        rdict[l]["theta"], rdict[l]["Qpol"], "o", color="C%i" % ccount
    )
    ax[1][0].plot(
        rdict[l]["theta"], rdict[l]["Qpolref"], "-", color="C%i" % ccount
    )
    ax[1][1].semilogy(
        rdict[l]["theta"],
        reldiff(rdict[l]["Qpol"], rdict[l]["Qpolref"]),
        "-",
        color="C%i" % ccount,
    )
    ax[2][0].plot(
        rdict[l]["theta"], rdict[l]["Qcpol"], "o", color="C%i" % ccount
    )
    ax[2][0].plot(
        rdict[l]["theta"], rdict[l]["Qcpolref"], "-", color="C%i" % ccount
    )
    ax[2][1].semilogy(
        rdict[l]["theta"],
        reldiff(rdict[l]["Qcpol"], rdict[l]["Qcpolref"]),
        "-",
        color="C%i" % ccount,
        label="$\\lambda{{}} = {0:.2f} \\mu{{}}$m".format(l),
    )
    ccount += 1

ax[2][1].legend(loc="best")
ax[0][0].set_ylabel("$Q_{ext}$")
ax[1][0].set_ylabel("$Q_{pol}$")
ax[2][0].set_ylabel("$Q_{cpol}$")
ax[0][1].set_ylabel("reldiff")
ax[1][1].set_ylabel("reldiff")
ax[2][1].set_ylabel("reldiff")
ax[2][0].set_xlabel("$\\theta{}$")
ax[2][1].set_xlabel("$\\theta{}$")

pl.tight_layout()
pl.savefig("1991Mishchenko.png", dpi=300, bbox_inches="tight")
