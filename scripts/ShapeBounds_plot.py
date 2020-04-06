import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 3.0)

thetas = np.linspace(0.0, 0.5 * np.pi, 50)

samples = np.arange(0.9, 0.98, 0.01)
cm = pl.get_cmap("viridis")
cb = pl.contourf(
    [[0.0, 0.0], [0.0, 0.0]], np.linspace(samples[0], samples[-1], 500), cmap=cm
)
pl.clf()

colwidth = 5
ax = [
    pl.subplot2grid((2, colwidth + 1), (0, 0), colspan=colwidth),
    pl.subplot2grid((2, colwidth + 1), (1, 0), colspan=colwidth),
    pl.subplot2grid((2, colwidth + 1), (0, colwidth), rowspan=2),
]
ref = None
ns = []
reldiffsQabs = []
reldiffsQabspol = []
maxdiffsQabs = []
maxdiffsQabspol = []
for isample in range(len(samples)):
    isample = len(samples) - isample - 1
    output = np.fromfile("ShapeBounds_i{0:02}.dat".format(isample)).reshape(
        (len(thetas), 2)
    )

    if ref is None:
        ref = output
    else:
        ns.append(samples[isample])
        reldiffsQabs.append(np.abs(output[:, 0] - ref[:, 0]).mean())
        reldiffsQabspol.append(np.abs(output[:, 1] - ref[:, 1]).mean())
        maxdiffsQabs.append(np.abs(output[:, 0] - ref[:, 0]).max())
        maxdiffsQabspol.append(np.abs(output[:, 1] - ref[:, 1]).max())

    ax[0].plot(
        thetas,
        output[:, 0],
        color=cm((samples[isample] - samples[0]) / (samples[-1] - samples[0])),
    )
    ax[1].plot(
        thetas,
        output[:, 1],
        color=cm((samples[isample] - samples[0]) / (samples[-1] - samples[0])),
    )

cbar = pl.colorbar(cb, cax=ax[2], label="sample fraction")
cbar.set_ticks(np.arange(0.9, 0.98, 0.01))

ax[0].set_xticklabels([])
ax[1].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[1].set_xticklabels(["$0$", "$\pi{}/4$", "$\pi{}/2$"])
ax[1].set_xlabel("$\\theta{}$")
ax[0].set_ylabel("$Q_{abs}$")
ax[1].set_ylabel("$Q_{abs,pol}$")
pl.tight_layout()
pl.savefig("ShapeBounds_fraction.png", dpi=300, bbox_inches="tight")
pl.savefig("ShapeBounds_fraction.pdf", dpi=300, bbox_inches="tight")
pl.close()

pl.rcParams["figure.figsize"] = (4.5, 2.5)

pl.semilogy(
    ns, reldiffsQabs / ref[:, 0].mean(), "-", color="C0", label="$Q_{abs}$"
)
pl.semilogy(
    ns,
    reldiffsQabspol / abs(ref[:, 1].mean()),
    "-",
    color="C1",
    label="$Q_{abs,pol}$",
)
pl.semilogy(ns, maxdiffsQabs / ref[:, 0].mean(), "--", color="C0")
pl.semilogy(ns, maxdiffsQabspol / abs(ref[:, 1].mean()), "--", color="C1")

meanline, = pl.plot([], [], "k-")
maxline, = pl.plot([], [], "k--")
leg = pl.legend([meanline, maxline], ["mean", "max"], loc="center")

pl.xlabel("sample fraction")
pl.ylabel("relative error")
pl.legend(loc="best")
pl.gca().add_artist(leg)
pl.tight_layout()
pl.savefig("ShapeBounds_fraction_convergence.png", dpi=300, bbox_inches="tight")
pl.savefig("ShapeBounds_fraction_convergence.pdf", dpi=300, bbox_inches="tight")
