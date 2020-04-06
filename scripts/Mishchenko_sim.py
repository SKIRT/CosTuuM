import numpy as np
import CosTuuM
import scipy.interpolate as interpol

mrdata = np.loadtxt("1985Draine.txt")
mr_lambda = mrdata[:, 0]
mr_nnorm = np.sqrt(mrdata[:, 1] ** 2 + mrdata[:, 2] ** 2)
mr_n = np.sqrt(0.5 * (mr_nnorm + mrdata[:, 1]))
mr_k = np.sqrt(0.5 * (mr_nnorm - mrdata[:, 1]))
mr_mr = mr_n + 1.0j * mr_k
mrspline = interpol.interp1d(mr_lambda, mr_mr, fill_value="extrapolate")


def get_refractive_index(wavelength, size, gtype):
    result = mrspline(wavelength * 1.0e6)
    return result


dust_properties = CosTuuM.CustomDustProperties(get_refractive_index)

nice_names = {
    "pdg": "Perfect Davis-Greenstein alignment",
    "ia": "Imperfect Mishchenko alignment",
    "ro": "Random orientations (no alignment)",
}
axis_ratios = {"prolate": 0.5, "oblate": 2.0}
orientations = {
    "pdg": {"prolate": 0.0, "oblate": 1.0},
    "ia": {"prolate": 1.0 / 5.0, "oblate": 3.0 / 5.0},
    "ro": {"prolate": 1.0 / 3.0, "oblate": 1.0 / 3.0},
}
r_ev = 2.0e-7
finethetas = np.linspace(0.0, np.pi / 2.0, 100)
wavelengths = np.array([0.2, 0.55, 1.0, 5.0]) * 1.0e-6


def test_case(axis_ratio_key, orientation_key):
    # we cheat to reproduce the different alignment behaviours:
    # for random orientations we set up an alignment distribution with a
    # transition size that is larger than the grain size
    if orientation_key == "ro":
        ad = CosTuuM.SizeBasedAlignmentDistribution(r_ev + 10.0, 0)
    elif orientation_key == "ia":
        ad = CosTuuM.SizeBasedAlignmentDistribution(
            r_ev - 10.0, CosTuuM.MISHCHENKO_ALIGNMENT
        )
    else:
        ad = CosTuuM.SizeBasedAlignmentDistribution(
            r_ev - 10.0, CosTuuM.DAVIS_GREENSTEIN_ALIGNMENT
        )
    sd = CosTuuM.SingleShapeShapeDistribution(axis_ratios[axis_ratio_key])

    results = CosTuuM.get_table(
        CosTuuM.SILICON,
        r_ev,
        wavelengths,
        finethetas,
        sd,
        ad,
        dust_properties,
        do_absorption=False,
        do_extinction=True,
    )
    results.tofile(
        "Mishchenko_{0}_{1}.dat".format(axis_ratio_key, orientation_key)
    )


for iar in range(len(axis_ratios.keys())):
    axis_ratio_key = [*axis_ratios][iar]
    for ior in range(len(orientations.keys())):
        orientation_key = [*orientations][ior]
        test_case(axis_ratio_key, orientation_key)
