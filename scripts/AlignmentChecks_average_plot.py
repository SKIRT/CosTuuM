import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (9.0, 3.5)

fig, ax = pl.subplots(2, 2, sharex=True, sharey="row")

wavelengths = np.logspace(-5.0, -3.0, 100) * 1.0e6
ds = [0.5, 2.0]
cos2beta = np.array([0.2, 0.6])
Rs = 2.0 * (1.5 * cos2beta - 0.5)
for id in range(len(ds)):
    d = ds[id]
    R = Rs[id]
    Qraw = np.fromfile("AlignmentChecks_d{0}_raw.dat".format(d)).reshape(
        (len(wavelengths), 2, 2)
    )
    Qpar = Qraw[:, 0, 0]
    Qperp = Qraw[:, 1, 0]

    Qabs = (2.0 * Qperp + Qpar) / 3.0
    Qabspol = -R * (Qpar - Qperp) / 3.0

    costheta, weight = np.polynomial.legendre.leggauss(50)
    theta = np.arccos(costheta)
    Qmish = np.fromfile("AlignmentChecks_d{0}_mish.dat".format(d)).reshape(
        (len(wavelengths), len(theta), 2)
    )
    Qmishquad = 0.5 * (weight[None, :, None] * Qmish).sum(axis=1)

    Qrand = np.fromfile("AlignmentChecks_d{0}_rand.dat".format(d)).reshape(
        (len(wavelengths), len(theta), 2)
    )
    Qperf = np.fromfile("AlignmentChecks_d{0}_perf.dat".format(d)).reshape(
        (len(wavelengths), len(theta), 2)
    )
    Qperfquad = 0.5 * (weight[None, :, None] * Qperf).sum(axis=1)

    label = None
    if id == 0:
        label = "picket fence"
    ax[0][id].loglog(wavelengths, Qabs, "-", label=label)
    ax[1][id].semilogx(wavelengths, np.abs(Qabspol) / Qabs)
    for f in [0.2, 0.4, 0.6]:
        Qmix = (1.0 - f) * Qrand[:, 0, :] + f * Qperfquad
        label = None
        if id == 1:
            label = "$f_A={0}$".format(f)
        ax[0][id].loglog(wavelengths, Qmix[:, 0], "-", label=label)
        ax[1][id].semilogx(wavelengths, np.abs(Qmix[:, 1]) / Qmix[:, 0], "-")

    label = None
    if id == 0:
        label = "CosTuuM"
    t = 3
    ax[0][id].loglog(wavelengths[::t], Qmishquad[::t, 0], ".", label=label)
    ax[1][id].semilogx(
        wavelengths[::t], np.abs(Qmishquad[::t, 1]) / Qmishquad[::t, 0], "."
    )

    ax[0][id].set_title("$d={0}$".format(d))

ax[0][0].set_ylim(5.0e-7, None)
ax[0][0].set_ylabel("$\\langle{}Q_{\\rm{}abs}\\rangle{}$")
ax[1][0].set_ylabel("$\\langle{}P_L\\rangle{}$")
ax[1][0].set_xlabel("$\\lambda{}$ ($\\mu{}$m)")
ax[1][1].set_xlabel("$\\lambda{}$ ($\\mu{}$m)")
ax[0][0].legend(loc="lower left")
ax[0][1].legend(
    loc="lower left", title="Dyck \& Beichman", ncol=3, columnspacing=1.0
)
pl.tight_layout()
pl.savefig("AlignmentChecks_average.png", dpi=300, bbox_inches="tight")
pl.savefig("AlignmentChecks_average.pdf", dpi=300, bbox_inches="tight")
