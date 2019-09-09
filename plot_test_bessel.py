import numpy as np
import scipy.special as special
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

testdata = np.loadtxt("test_bessel.txt")

refx = np.linspace(0.05, 100.0, 2000)
refj = special.spherical_jn(10, refx + 1.0j * refx)
refdj = special.spherical_jn(10, refx + 1.0j * refx, derivative=True)

fig, ax = pl.subplots(2, 2)

ax[0][0].plot(testdata[:, 0], testdata[:, 1], ".")
ax[0][0].semilogy(refx, refj.real)
ax[0][1].plot(testdata[:, 0], testdata[:, 2], ".")
ax[0][1].semilogy(refx, refj.imag)

ax[1][0].plot(testdata[:, 0], testdata[:, 3], ".")
ax[1][0].semilogy(refx, refdj.real)
ax[1][1].plot(testdata[:, 0], testdata[:, 4], ".")
ax[1][1].semilogy(refx, refdj.imag)

pl.tight_layout()
pl.savefig("test_bessel.png", dpi=300)
