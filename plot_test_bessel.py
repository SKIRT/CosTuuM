##
# @file plot_test_bessel.py
#
# @brief Script to plot the output from the testBesselFunctions unit test.
#
# @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
##

import numpy as np
import scipy.special as special
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

# maximum order (should match the value in the unit test)
nmax = 80

# open the complex output from the unit test
testdata = np.loadtxt("test_bessel_complex.txt")
# filter out the 1st and nth order, real and imaginary parts, and function and
# derivative
testj1r = testdata[::2, 0]
testj1i = testdata[::2, 1]
testjnr = testdata[1::2, 0]
testjni = testdata[1::2, 1]
testdj1r = testdata[::2, 2]
testdj1i = testdata[::2, 3]
testdjnr = testdata[1::2, 2]
testdjni = testdata[1::2, 3]

# set up the reference values
# the refx array should match the x expression in the unit test
refx = np.arange(0.05, 100.0, 0.1)
# compute references for the 1st and nth order complex test
refj1 = special.spherical_jn(1, refx + 1.0j * refx)
refdj1 = special.spherical_jn(1, refx + 1.0j * refx, derivative=True)
# convert the derivatives into the appropriate derivative expression
refdj1 += refj1 / (refx + 1.0j * refx)
refjn = special.spherical_jn(nmax, refx + 1.0j * refx)
refdjn = special.spherical_jn(nmax, refx + 1.0j * refx, derivative=True)
refdjn += refjn / (refx + 1.0j * refx)

# plot the complex results
fig, ax = pl.subplots(2, 2)

# functions themselves
ax[0][0].plot(refx, testj1r, ".")
ax[0][0].semilogy(refx, refj1.real)
ax[0][0].plot(refx, testjnr, ".")
ax[0][0].semilogy(refx, refjn.real)
ax[0][1].plot(refx, testj1i, ".")
ax[0][1].semilogy(refx, refj1.imag)
ax[0][1].plot(refx, testjni, ".")
ax[0][1].semilogy(refx, refjn.imag)

# derivatives
ax[0][0].plot(refx, testdj1r, ".")
ax[0][0].semilogy(refx, refdj1.real)
ax[0][0].plot(refx, testdjnr, ".")
ax[0][0].semilogy(refx, refdjn.real)
ax[0][1].plot(refx, testdj1i, ".")
ax[0][1].semilogy(refx, refdj1.imag)
ax[0][1].plot(refx, testdjni, ".")
ax[0][1].semilogy(refx, refdjn.imag)

# relative differences between test output and reference values
ax[1][0].semilogy(refx, abs(refj1.real - testj1r) / abs(refj1.real + testj1r))
ax[1][0].semilogy(refx, abs(refjn.real - testjnr) / abs(refjn.real + testjnr))
ax[1][0].semilogy(
    refx, abs(refdj1.real - testdj1r) / abs(refdj1.real + testdj1r)
)
ax[1][0].semilogy(
    refx, abs(refdjn.real - testdjnr) / abs(refdjn.real + testdjnr)
)
ax[1][1].semilogy(refx, abs(refj1.imag - testj1i) / abs(refj1.imag + testj1i))
ax[1][1].semilogy(refx, abs(refjn.imag - testjni) / abs(refjn.imag + testjni))
ax[1][1].semilogy(
    refx, abs(refdj1.imag - testdj1i) / abs(refdj1.imag + testdj1i)
)
ax[1][1].semilogy(
    refx, abs(refdjn.imag - testdjni) / abs(refdjn.imag + testdjni)
)

# save complex figure and close plot
pl.tight_layout()
pl.savefig("test_bessel_complex.png", dpi=300)
pl.close()

# repeat for the real test...
testdata = np.loadtxt("test_bessel_real.txt")
# filter out 1st and nth order, 1st and 2nd kind, functions and derivatives
testj1 = testdata[::2, 0]
testy1 = testdata[::2, 2]
testjn = testdata[1::2, 0]
testyn = testdata[1::2, 2]
testdj1 = testdata[::2, 1]
testdy1 = testdata[::2, 3]
testdjn = testdata[1::2, 1]
testdyn = testdata[1::2, 3]

# set up reference values
refj1 = special.spherical_jn(1, refx)
refdj1 = special.spherical_jn(1, refx, derivative=True)
# apply correction for derivative expression
refdj1 += refj1 / refx
refjn = special.spherical_jn(nmax, refx)
refdjn = special.spherical_jn(nmax, refx, derivative=True)
refdjn += refjn / refx
refy1 = special.spherical_yn(1, refx)
refdy1 = special.spherical_yn(1, refx, derivative=True)
refdy1 += refy1 / refx
refyn = special.spherical_yn(nmax, refx)
refdyn = special.spherical_yn(nmax, refx, derivative=True)
refdyn += refyn / refx

# plot results
fig, ax = pl.subplots(2, 2)

# functions
ax[0][0].plot(refx, testj1, ".")
ax[0][0].semilogy(refx, refj1)
ax[0][0].plot(refx, testjn, ".")
ax[0][0].semilogy(refx, refjn)
ax[0][1].plot(refx, testy1, ".")
ax[0][1].semilogy(refx, refy1)
ax[0][1].plot(refx, testyn, ".")
ax[0][1].semilogy(refx, refyn)

# derivatives
ax[0][0].plot(refx, testdj1, ".")
ax[0][0].semilogy(refx, refdj1)
ax[0][0].plot(refx, testdjn, ".")
ax[0][0].semilogy(refx, refdjn)
ax[0][1].plot(refx, testdy1, ".")
ax[0][1].semilogy(refx, refdy1)
ax[0][1].plot(refx, testdyn, ".")
ax[0][1].semilogy(refx, refdyn)

# relative differences between test output and reference values
ax[1][0].semilogy(refx, abs(refj1 - testj1) / abs(refj1 + testj1))
ax[1][0].semilogy(refx, abs(refjn - testjn) / abs(refjn + testjn))
ax[1][0].semilogy(refx, abs(refdj1 - testdj1) / abs(refdj1 + testdj1))
ax[1][0].semilogy(refx, abs(refdjn - testdjn) / abs(refdjn + testdjn))
ax[1][1].semilogy(refx, abs(refy1 - testy1) / abs(refy1 + testy1))
ax[1][1].semilogy(refx, abs(refyn - testyn) / abs(refyn + testyn))
ax[1][1].semilogy(refx, abs(refdy1 - testdy1) / abs(refdy1 + testdy1))
ax[1][1].semilogy(refx, abs(refdyn - testdyn) / abs(refdyn + testdyn))

# save real output plot
pl.tight_layout()
pl.savefig("test_bessel_real.png", dpi=300)
