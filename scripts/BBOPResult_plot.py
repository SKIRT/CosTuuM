import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 5.0)

theta = np.linspace(0.0, 0.5 * np.pi, 50)
wav = np.array([7.0e-5, 2.0e-4, 3.5e-4])


quad = np.fromfile("BBOPResult.dat").reshape((len(wav), len(theta), 2))

fig, ax = pl.subplots(3, 1, sharex=True)
ax[0].semilogy(theta, quad[:, :, 0].T)
ax[1].semilogy(theta[1:], -quad[:, 1:, 1].T)
ax[2].plot(theta, np.abs(quad[0, :, 1]) / quad[0, :, 0], label="$70~\\mu{}$m")
ax[2].plot(theta, np.abs(quad[1, :, 1]) / quad[1, :, 0], label="$200~\\mu{}$m")
ax[2].plot(theta, np.abs(quad[2, :, 1]) / quad[2, :, 0], label="$350~\\mu{}$m")

ax[0].set_ylabel("$\\langle{}K_{a,I}\\rangle{}$ (code units)")
ax[1].set_ylabel("$-\\langle{}K_{a,Q}\\rangle{}$ (code units)")
ax[2].set_ylabel("$P_L$")
ax[2].set_xlabel("$\\theta{}$")
ax[2].set_xticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
ax[2].set_xticklabels(["$0$", "$\\pi{}/4$", "$\\pi{}/2$"])

ax[2].legend(loc="best")
pl.tight_layout()
pl.savefig("FullPlot.png", dpi=300, bbox_inches="tight")
pl.savefig("FullPlot.pdf", dpi=300, bbox_inches="tight")
