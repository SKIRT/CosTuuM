import numpy as np
import CosTuuM

ad = CosTuuM.SizeBasedAlignmentDistribution(1.0e-7, 0)
dp = CosTuuM.DraineDustProperties()
thetas = np.linspace(0.0, 0.5 * np.pi, 50)

nmin = 1
nmax = 20
for ngauss in range(nmin, nmax + 1):
    sd = CosTuuM.DraineHensleyShapeDistribution(ngauss, 0.0, 0.96)
    output = CosTuuM.get_table(
        types=CosTuuM.SILICON,
        sizes=1.0e-6,
        wavelengths=1.0e-5,
        thetas=thetas,
        shape_distribution=sd,
        alignment_distribution=ad,
        dust_properties=dp,
        minimum_order=10,
        maximum_order=146,
        gauss_legendre_factor=2,
        tolerance=1.0e-4,
        number_of_quadrature_angles=20,
        number_of_threads=1,
        verbose=True,
        account_for_scattering=True,
        maximum_memory_size=5000000000,
    )

    output.tofile("ShapeQuadrature_n{0:02d}.dat".format(ngauss))
