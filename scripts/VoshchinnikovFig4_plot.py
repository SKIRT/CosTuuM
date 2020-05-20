import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (9.0, 3.0)

vosh = np.load("1993Voshchinnikov_fig4.npz", mmap_mode="r")

mrs = vosh["mr"]
ds = vosh["d"]
thetas = vosh["theta"]
Qvosh = vosh["Qpol"]
rV = 1.5e-7
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

cost = np.load("VoshchinnikovFig4.npz", mmap_mode="r")
Qcos = cost["Qcos"]

colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7"]
fig, ax = pl.subplots(1, 3, sharex=True)
ax[0].set_title("prolate, $m_r=1.31+0.01i$")
for id in range(4):
    ax[0].plot(thetas, Qcos[0, id, :], "-", color=colors[id], alpha=0.5)
    ax[0].plot(thetas[::8], Qvosh[0, id, ::8], "x", color=colors[id])
ax[1].set_title("oblate, $m_r=1.31+0.01i$")
for id in range(4, 8):
    ax[1].plot(thetas, np.abs(Qcos[0, id, :]), "-", color=colors[id], alpha=0.5)
    ax[1].plot(thetas[::8], np.abs(Qvosh[0, id, ::8]), "x", color=colors[id])
ax[2].set_title("prolate, $m_r=1.70+0.03i$")
for id in range(4):
    ax[2].plot(thetas, Qcos[1, id, :], "-", color=colors[id], alpha=0.5)
    ax[2].plot(thetas[::8], Qvosh[1, id, ::8], "x", color=colors[id])

ax[0].set_xticks([0, 0.25 * np.pi, 0.5 * np.pi])
ax[0].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])

ax[0].set_ylabel("$|Q_{ext,pol}|$")

ax[1].set_ylim(-0.01, 0.3)
ax[0].set_ylim(ax[1].get_ylim())
ax[1].set_yticklabels([])

ax[2].set_ylim(-0.01, 1.0)

ax[0].set_xlabel("$\\theta{}$")
ax[1].set_xlabel("$\\theta{}$")
ax[2].set_xlabel("$\\theta{}$")

for id in range(4):
    ax[2].plot(
        [], [], "-", color=colors[id], label="$d={0:.2f}$".format(ds[id])
    )
ax[2].legend(loc="best", handlelength=0.5, labelspacing=0.1)
for id in range(4, 8):
    ax[1].plot(
        [], [], "-", color=colors[id], label="$d={0:.1f}$".format(ds[id])
    )
ax[1].legend(loc="best", ncol=2, handlelength=0.5, labelspacing=0.1)

ax[0].plot([], [], "k-", label="\\textsc{CosTuuM}")
ax[0].plot([], [], "kx", label="Voshchinnikov \& Faranov (1993)")
ax[0].legend(loc="best", handlelength=0.5, labelspacing=0.1)

pl.tight_layout()
pl.savefig("VoshchinnikovFig4.pdf", dpi=300, bbox_inches="tight")
