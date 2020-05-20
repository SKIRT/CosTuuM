import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (9.0, 3.0)

vosh = np.load("1993Voshchinnikov_fig3.npz", mmap_mode="r")

mrs = vosh["mr"]
ds = vosh["d"]
rVs = vosh["rV"]
thetas = vosh["theta"]
Qvosh = vosh["Qext"]
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

cost = np.load("VoshchinnikovFig3.npz", mmap_mode="r")
Qcos = cost["Qcos"]


colors = ["C0", "C1", "C2", "C3", "C4"]
fig, ax = pl.subplots(1, 3, sharex=True)
ax[0].set_title("$d=0.5, m_r=1.31+0.01i$")
for irV in range(len(rVs)):
    ax[0].plot(thetas, Qcos[0, 0, irV, :], "-", color=colors[irV], alpha=0.5)
    ax[0].plot(thetas[::8], Qvosh[0, 0, irV, ::8], "x", color=colors[irV])
ax[1].set_title("$d=2.0, m_r=1.31+0.01i$")
for irV in range(len(rVs)):
    ax[1].plot(thetas, Qcos[0, 1, irV, :], "-", color=colors[irV], alpha=0.5)
    ax[1].plot(thetas[::8], Qvosh[0, 1, irV, ::8], "x", color=colors[irV])
ax[2].set_title("$d=0.5, m_r=1.70+0.03i$")
for irV in range(len(rVs)):
    ax[2].plot(thetas, Qcos[1, 0, irV, :], "-", color=colors[irV], alpha=0.5)
    ax[2].plot(thetas[::8], Qvosh[1, 0, irV, ::8], "x", color=colors[irV])

ax[0].set_xticks([0, 0.25 * np.pi, 0.5 * np.pi])
ax[0].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])

ax[0].set_ylabel("$Q_{ext}$")

ax[0].set_xlabel("$\\theta{}$")
ax[1].set_xlabel("$\\theta{}$")
ax[2].set_xlabel("$\\theta{}$")

for irV in range(len(rVs)):
    ax[1].plot(
        [],
        [],
        "-",
        color=colors[irV],
        label="$a={0:.2f}~\\mu{{}}$m".format(rVs[irV] * 1.0e6),
    )
ax[1].legend(loc="best", ncol=2, handlelength=0.5, labelspacing=0.1)

ax[0].plot([], [], "k-", label="\\textsc{CosTuuM}")
ax[0].plot([], [], "kx", label="Voshchinnikov \& Farafonov (1993)")
ax[0].legend(loc="best", handlelength=0.5, labelspacing=0.1)

ax[0].set_ylim(-0.1, 3.5)
ax[1].set_ylim(-0.1, 3.5)
ax[1].set_yticklabels([])

pl.tight_layout()
pl.savefig("VoshchinnikovFig3.pdf", dpi=300, bbox_inches="tight")
