/**
 * @file CTM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Matrix.hpp"
#include "SpecialFunctions.hpp"

#include <cinttypes>
#include <cmath>
#include <complex>
#include <vector>

/**
 * @brief Main program entry point.
 *
 * For now, does the entire T-matrix calculation for the Mishchenko default
 * parameters.
 *
 * @param argc Number of command line arguments (currently ignored).
 * @param argv Command line arguments (currently ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// input parameters (should become real parameters at some point
  // size of the particle (in same units as the wavelength)
  const double axi = 10.;
  // ratio between the equal surface area sphere radius and equal volume sphere
  // radius (is recomputed if not equal to 1)
  double ratio_of_radii = 0.1;
  // wavelength of incoming radiation (in same units as the particle size)
  const double wavelength = 2. * M_PI;
  // refractive index
  const std::complex<double> mr(1.5, 0.02);
  // ratio of horizontal and vertical axis of the spheroidal particle
  const double axis_ratio = 0.5;
  // tolerance for the calculation
  const double tolerance = 1.e-4;
  // number of Gauss-Legendre points to use as a multiplier of the maximum
  // order of spherical harmonics used to decompose the electric field, during
  // the first loop of the algorithm
  const uint_fast32_t ndgs = 1;

  /// hardcoded program parameters
  // maximum number of iterations during the first loop
  const uint_fast32_t maximum_order = 200;
  // maximum number of iterations during the second loop
  const uint_fast32_t npng1 = 500;

  (void)npng1;

  // make sure 'rat' contains the right ratio if it is not 1
  if (std::abs(ratio_of_radii - 1.) > 1.e-8) {
    ratio_of_radii =
        SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
            axis_ratio);
  }

  // R_V is the equivalent sphere radius
  const double R_V = ratio_of_radii * axi;
  // we need a reasonable initial guess for the order of the spherical harmonics
  // the below expression provides this guess
  const double xev = 2. * M_PI * R_V / wavelength;
  uint_fast32_t nmax = std::max(4., xev + 4.05 * std::cbrt(xev));

  // loop control variables
  double old_qext = 0.;
  double old_qsca = 0.;
  double dext = 1.;
  double dsca = 1.;
  while (nmax < maximum_order && (dext > tolerance || dsca > tolerance)) {
    const uint_fast32_t ngauss = ndgs * nmax;
    const uint_fast32_t ngauss2 = 2 * ngauss;

    // this is basically func_const() in the Python version
    std::vector<double> an(nmax);
    Matrix<double> ann(nmax, nmax);
    std::vector<double> dd(nmax);
    std::vector<double> sinthetainv(ngauss2);
    std::vector<double> sintheta2inv(ngauss2);
    std::vector<double> costheta(ngauss2), wgauss(ngauss2);
    for (uint_fast32_t ni = 0; ni < nmax; ++ni) {
      const double nn = (ni + 2.) * (ni + 1.);
      an[ni] = nn;
      const double d = std::sqrt((2. * (ni + 1.) + 1.) / nn);
      dd[ni] = d;
      for (uint_fast32_t nj = 0; nj < ni + 1; ++nj) {
        const double ddd = 0.5 * d * dd[ni];
        ann(ni, nj) = ddd;
        ann(nj, ni) = ddd;
      }
    }
    SpecialFunctions::get_gauss_legendre_points_and_weigths(ngauss2, costheta,
                                                            wgauss);
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      const double this_sintheta2inv = 1. / (1. - costheta[ig] * costheta[ig]);
      sintheta2inv[ig] = this_sintheta2inv;
      sintheta2inv[ngauss2 - ig - 1] = this_sintheta2inv;
      const double this_sinthetainv = std::sqrt(this_sintheta2inv);
      sinthetainv[ig] = this_sinthetainv;
      sinthetainv[ngauss2 - ig - 1] = this_sinthetainv;
    }

    // this is func_vary() in the Python version
    std::vector<double> krinv(ngauss2);
    std::vector<std::complex<double>> krmrinv(ngauss2);
    std::vector<double> kr(ngauss2);
    std::vector<std::complex<double>> krmr(ngauss2);
    std::vector<double> r2(ngauss2);
    std::vector<double> dr_over_r(ngauss2);
    SpecialFunctions::get_r_dr_spheroid(costheta, R_V, axis_ratio, r2,
                                        dr_over_r);
    const double wavenumber = 2. * M_PI / wavelength;
    const double wavenumber2 = wavenumber * wavenumber;
    const std::complex<double> wavenumber2_mr = mr * wavenumber;
    const std::complex<double> mrinv = 1. / mr;
    for (uint_fast32_t i = 0; i < ngauss2; ++i) {
      const double r = std::sqrt(r2[i]);
      kr[i] = wavenumber * r;
      krmr[i] = mr * kr[i];
      krinv[i] = 1. / kr[i];
      krmrinv[i] = mrinv * krinv[i];
    }
    Matrix<double> jkr(ngauss2, nmax), djkr(ngauss2, nmax), ykr(ngauss2, nmax),
        dykr(ngauss2, nmax);
    Matrix<std::complex<double>> jkrmr(ngauss2, nmax), djkrmr(ngauss2, nmax);
    for (uint_fast32_t ig = 0; ig < ngauss2; ++ig) {
      SpecialFunctions::spherical_j_jdj_array(nmax, kr[ig], jkr.get_row(ig),
                                              djkr.get_row(ig));
      SpecialFunctions::spherical_y_ydy_array(nmax, kr[ig], ykr.get_row(ig),
                                              dykr.get_row(ig));
      SpecialFunctions::spherical_j_jdj_array(nmax, krmr[ig], jkrmr.get_row(ig),
                                              djkrmr.get_row(ig));
    }

    // func_tmatr0()
    const uint_fast32_t nmax2 = 2 * nmax;

    std::vector<int_fast8_t> signs(nmax2);
    Matrix<double> werner_d(ngauss2, nmax);
    Matrix<double> dwerner_d(ngauss2, nmax);
    std::vector<double> wr2(ngauss);
    Matrix<std::complex<double>> J12(nmax, nmax);
    Matrix<std::complex<double>> J21(nmax, nmax);
    Matrix<std::complex<double>> RgJ12(nmax, nmax);
    Matrix<std::complex<double>> RgJ21(nmax, nmax);
    Matrix<std::complex<double>> Q(nmax2, nmax2);
    Matrix<std::complex<double>> RgQ(nmax2, nmax2);

    int_fast8_t si = 1;
    for (uint_fast32_t m = 0; m < nmax2; ++m) {
      si = -si;
      signs[m] = si;
    }
    for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
      const uint_fast32_t i1 = ngauss + ig;
      const uint_fast32_t i2 = ngauss - ig + 1;
      std::vector<double> dv1(nmax), dv2(nmax);
      SpecialFunctions::wigner_dn_0m(costheta[ig], nmax, 0, &dv1[0], &dv2[0]);
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        si = signs[n];
        werner_d(i1 - 1, n) = dv1[n];
        werner_d(i2 - 1, n) = si * dv1[n];
        dwerner_d(i1 - 1, n) = dv2[n];
        dwerner_d(i2 - 1, n) = -si * dv2[n];
      }
    }
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      wr2[ig] = wgauss[ig] * r2[ig];
    }
    // compute J and RgJ components...

    (void)wavenumber2;
    (void)wavenumber2_mr;
    (void)old_qext;
    (void)old_qsca;
    ++nmax;
  }

  return 0;
}
