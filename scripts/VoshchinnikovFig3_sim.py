import numpy as np
import CosTuuM

vosh = np.load("1993Voshchinnikov_fig3.npz", mmap_mode="r")

mrs = vosh["mr"]
ds = vosh["d"]
rVs = vosh["rV"]
thetas = vosh["theta"]
Qvosh = vosh["Qext"]
wav = 5.0e-7

# convert theta to radians
thetas *= np.pi / 180.0

Qcos = np.zeros((mrs.shape[0], ds.shape[0], rVs.shape[0], thetas.shape[0]))
for imr in range(len(mrs)):
    mr = mrs[imr]
    print("mr:", mr)
    for id in range(len(ds)):
        d = ds[id]
        print("d:", d)

        def get_refractive_index(wavelength, size, gtype):
            return np.array(mr)

        dust_properties = CosTuuM.CustomDustProperties(get_refractive_index)
        results = CosTuuM.get_table(
            CosTuuM.SILICON,
            rVs,
            wav,
            thetas,
            CosTuuM.SingleShapeShapeDistribution(d),
            CosTuuM.SizeBasedAlignmentDistribution(
                rVs[0] - 10.0, CosTuuM.DISABLE_ALIGNMENT
            ),
            dust_properties,
            do_extinction=True,
            do_absorption=False,
        )

        Qcos[imr, id, :, :] = results[:, :, 0]


np.savez("VoshchinnikovFig3.npz", theta=thetas, Qcos=Qcos)
