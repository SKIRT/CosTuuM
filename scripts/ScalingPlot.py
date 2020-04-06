import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (4.5, 3.5)

fig, ax = pl.subplots(2, 1, sharex=True, sharey=True)

names = [
    "timings_stormy_bind.dat",
    "timings_stormy_nobind.dat",
    "timings_gandalf_bind.dat",
    "timings_gandalf_nobind.dat",
]
colors = ["C0", "C0", "C1", "C1"]
lines = ["-", "--", "-", "--"]
nthreads = [56, 56, 64, 64]
for i in range(len(names)):
    threads = np.arange(1, nthreads[i] + 1)
    times = np.fromfile(names[i])
    speedup = times[0] / times
    fakespeedup = 2 * times[1] / times
    pe = speedup / threads
    fakepe = fakespeedup / threads

    ax[0].plot(threads, pe, lines[i], color=colors[i])
    ax[1].plot(threads[1:], fakepe[1:], lines[i], color=colors[i])

ax[0].plot([], [], color=colors[0], label="Xeon E5-2690")
ax[0].plot([], [], color=colors[2], label="Xeon E5-4650")
ax[1].plot([], [], "k-", label="pinning")
ax[1].plot([], [], "k--", label="no pinning")

ax[0].plot(threads, np.ones(threads.shape), ":", color="C2")
ax[1].plot(threads, np.ones(threads.shape), ":", color="C2")

ax[0].legend(loc="upper right")
ax[1].legend(loc="upper right")

ax[0].set_ylabel("$(1/n) t(1)/t(n)$")
ax[1].set_ylabel("$(2/n) t(2)/t(n)$")

ax[1].set_xlabel("number of threads")

ax[1].set_ylim(0.0, 1.2)

pl.tight_layout()
pl.savefig("ScalingTest.png", dpi=300, bbox_inches="tight")
pl.savefig("ScalingTest.pdf", dpi=300, bbox_inches="tight")
