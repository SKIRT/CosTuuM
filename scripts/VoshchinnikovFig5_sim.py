import numpy as np
import CosTuuM

vosh = np.load("1993Voshchinnikov_fig5.npz", mmap_mode="r")

mr = 1.31 + 0.01j
ds = vosh["d"]
thetas = vosh["theta"]
Qs = vosh["Qs"]
rV = 1.5e-7
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

Qs = np.zeros((ds.shape[0], thetas.shape[0], 2))
for id in range(len(ds)):
    d = ds[id]
    print("d:", d)

    def get_refractive_index(wavelength, size, gtype):
        return np.array(mr)

    dust_properties = CosTuuM.CustomDustProperties(get_refractive_index)
    results = CosTuuM.get_table(
        CosTuuM.SILICON,
        rV,
        wav,
        thetas,
        CosTuuM.SingleShapeShapeDistribution(d),
        CosTuuM.SizeBasedAlignmentDistribution(
            rV - 10.0, CosTuuM.DISABLE_ALIGNMENT
        ),
        dust_properties,
        do_extinction=True,
        do_absorption=False,
    )

    Qs[id, :, 0] = results[:, 0]
    Qs[id, :, 1] = results[:, 1]


np.savez("VoshchinnikovFig5.npz", theta=thetas, Qs=Qs)
