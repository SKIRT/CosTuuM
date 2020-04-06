import numpy as np
import CosTuuM

amin = 5.0e-9
amax = 2.0e-6

ad = CosTuuM.SizeBasedAlignmentDistribution(1.0e-7, 0)
dp = CosTuuM.DraineDustProperties()

theta = np.linspace(0.0, 0.5 * np.pi, 50)
wav = np.array([7.0e-5, 2.0e-4, 3.5e-4])


def integrand(a):
    output = CosTuuM.get_table(
        types=CosTuuM.SILICON,
        sizes=a,
        wavelengths=wav,
        thetas=theta,
        shape_distribution=CosTuuM.DraineHensleyShapeDistribution(
            20, 0.0, 0.96
        ),
        alignment_distribution=ad,
        dust_properties=dp,
        minimum_order=10,
        maximum_order=146,
        gauss_legendre_factor=2,
        tolerance=1.0e-4,
        number_of_quadrature_angles=20,
        number_of_threads=4,
        verbose=True,
        account_for_scattering=True,
        maximum_memory_size=30000000000,
    )
    return a[:, None, None, None] ** (-1.5) * np.pi * output


abreak = 1.0e-7
wfac_small = 0.5 * (abreak - amin)
xterm_small = 0.5 * (abreak + amin)
wfac_large = 0.5 * (amax - abreak)
xterm_large = 0.5 * (amax + abreak)
ngauss = 128
ag, wg = np.polynomial.legendre.leggauss(ngauss // 2)
ag_small = wfac_small * ag + xterm_small
wg_small = wfac_small * wg
ag_large = wfac_large * ag + xterm_large
wg_large = wfac_large * wg
intasmall = wg_small[:, None, None, None] * integrand(ag_small)
intalarge = wg_large[:, None, None, None] * integrand(ag_large)
quad = intasmall.sum(axis=0) + intalarge.sum(axis=0)

quad.tofile("BBOPResult.dat")
