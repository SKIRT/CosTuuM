################################################################################
# This file is part of CosTuuM
# Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CosTuuM is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
###############################################################################

##
# @file plot_test_special_functions.py
#
# @brief Script to plot the output from the testSpecialFunctions unit test.
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

### Bessel functions

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
pl.close()

### Wigner D function

# open the output from the unit test
testdata = np.loadtxt("test_wigner_d.txt")
# filter out different m values and functions and derivatives
testdm0 = testdata[::2, 0]
testddm0 = testdata[::2, 1]
testdmn = testdata[1::2, 0]
testddmn = testdata[1::2, 1]

# set up the reference range (should match the expression in the unit test)
refcosx = np.arange(-0.999, 1.0, 0.002)

# generate reference values for m=0
m = 0
refdm0 = np.zeros(refcosx.shape)
refddm0 = np.zeros(refcosx.shape)
# special.lpmn does not take input arrays for x, so we need to iterate
for i in range(len(refcosx)):
    lpmn, dlpmn = special.lpmn(m, nmax + 1, refcosx[i])
    refdm0[i] = lpmn[m, nmax]
    refddm0[i] = dlpmn[m, nmax]
# compute the conversion factor from P^m_n to d^n_{0m}
factor = (-1.0) ** m * np.sqrt(
    np.math.factorial(nmax - m) / np.math.factorial(nmax + m)
)
refdm0 *= factor
# also add the correction factor for the derivative:
# d/dx P(cos(x)) = d/dz P(z) d/dx cos(x) = -sin(x) d/dz P(z)
refddm0 *= -factor * np.sqrt(1.0 - refcosx ** 2)

# generate reference values for m=nmax
m = nmax
refdmn = np.zeros(refcosx.shape)
refddmn = np.zeros(refcosx.shape)
for i in range(len(refcosx)):
    lpmn, dlpmn = special.lpmn(m, nmax + 1, refcosx[i])
    refdmn[i] = lpmn[m, nmax]
    refddmn[i] = dlpmn[m, nmax]
factor = (-1.0) ** m * np.sqrt(
    np.math.factorial(nmax - m) / np.math.factorial(nmax + m)
)
refdmn *= factor
refddmn *= -factor * np.sqrt(1.0 - refcosx ** 2)

# plot the results
fig, ax = pl.subplots(2, 2)

# plot functions and derivatives
ax[0][0].plot(refcosx, testdm0, ".")
ax[0][0].plot(refcosx, refdm0)
ax[0][0].plot(refcosx, testdmn, ".")
ax[0][0].plot(refcosx, refdmn)
ax[0][1].plot(refcosx, testddm0, ".")
ax[0][1].plot(refcosx, refddm0)
ax[0][1].plot(refcosx, testddmn, ".")
ax[0][1].plot(refcosx, refddmn)

# plot the relative difference between the reference and the test results
ax[1][0].semilogy(refcosx, abs(refdm0 - testdm0) / abs(refdm0 + testdm0))
ax[1][0].semilogy(refcosx, abs(refdmn - testdmn) / abs(refdmn + testdmn))
ax[1][1].semilogy(refcosx, abs(refddm0 - testddm0) / abs(refddm0 + testddm0))
ax[1][1].semilogy(refcosx, abs(refddmn - testddmn) / abs(refddmn + testddmn))

# save the figure
pl.tight_layout()
pl.savefig("test_wigner_d.png", dpi=300)
pl.close()

# open the output from the unit test
testdata = np.loadtxt("test_wigner_d_sinx.txt")
# filter out different m values and functions and derivatives
testdm0 = testdata[::2, 0]
testddm0 = testdata[::2, 1]
testdmn = testdata[1::2, 0]
testddmn = testdata[1::2, 1]

# set up the reference range (should match the expression in the unit test)
refcosx = np.arange(-0.999, 1.0, 0.002)

# generate reference values for m=0
m = 0
refdm0 = np.zeros(refcosx.shape)
refddm0 = np.zeros(refcosx.shape)
# special.lpmn does not take input arrays for x, so we need to iterate
for i in range(len(refcosx)):
    lpmn, dlpmn = special.lpmn(m, nmax + 1, refcosx[i])
    refdm0[i] = lpmn[m, nmax]
    refddm0[i] = dlpmn[m, nmax]
# compute the conversion factor from P^m_n to d^n_{0m}
factor = (-1.0) ** m * np.sqrt(
    np.math.factorial(nmax - m) / np.math.factorial(nmax + m)
)
refdm0 *= factor
# also add the correction factor for the derivative:
# d/dx P(cos(x)) = d/dz P(z) d/dx cos(x) = -sin(x) d/dz P(z)
refddm0 *= -factor * np.sqrt(1.0 - refcosx ** 2)

refdm0 /= np.sqrt(1.0 - refcosx ** 2)

# generate reference values for m=nmax
m = nmax
refdmn = np.zeros(refcosx.shape)
refddmn = np.zeros(refcosx.shape)
for i in range(len(refcosx)):
    lpmn, dlpmn = special.lpmn(m, nmax + 1, refcosx[i])
    refdmn[i] = lpmn[m, nmax]
    refddmn[i] = dlpmn[m, nmax]
factor = (-1.0) ** m * np.sqrt(
    np.math.factorial(nmax - m) / np.math.factorial(nmax + m)
)
refdmn *= factor
refddmn *= -factor * np.sqrt(1.0 - refcosx ** 2)

refdmn /= np.sqrt(1.0 - refcosx ** 2)

# plot the results
fig, ax = pl.subplots(2, 2)

# plot functions and derivatives
ax[0][0].plot(refcosx, testdm0, ".")
ax[0][0].plot(refcosx, refdm0)
ax[0][0].plot(refcosx, testdmn, ".")
ax[0][0].plot(refcosx, refdmn)
ax[0][1].plot(refcosx, testddm0, ".")
ax[0][1].plot(refcosx, refddm0)
ax[0][1].plot(refcosx, testddmn, ".")
ax[0][1].plot(refcosx, refddmn)

# plot the relative difference between the reference and the test results
ax[1][0].semilogy(refcosx, abs(refdm0 - testdm0) / abs(refdm0 + testdm0))
ax[1][0].semilogy(refcosx, abs(refdmn - testdmn) / abs(refdmn + testdmn))
ax[1][1].semilogy(refcosx, abs(refddm0 - testddm0) / abs(refddm0 + testddm0))
ax[1][1].semilogy(refcosx, abs(refddmn - testddmn) / abs(refddmn + testddmn))

# save the figure
pl.tight_layout()
pl.savefig("test_wigner_d_sinx.png", dpi=300)
pl.close()

testdata = np.loadtxt("test_gauss_legendre_quadrature.txt")

refx, refw = np.polynomial.legendre.leggauss(200)

fig, ax = pl.subplots(2, 1)

ax[0].plot(testdata[:, 0], testdata[:, -1], ".")
ax[0].plot(refx, refw)

ax[1].semilogy(refx, abs(refx - testdata[:, 0]) / abs(refx + testdata[:, 0]))
ax[1].semilogy(refx, abs(refw - testdata[:, 1]) / abs(refw + testdata[:, 1]))

pl.tight_layout()
pl.savefig("test_gauss_legendre_quadrature.png", dpi=300)
