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
    Matrix<double> ann(nmax);
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
    std::vector<double> jkr(ngauss2), djkr(ngauss2), ykr(ngauss2),
        dykr(ngauss2);
    std::vector<std::complex<double>> jkrmr(ngauss2), djkrmr(ngauss2);
    for (uint_fast32_t ig = 0; ig < ngauss2; ++ig) {
      //        SpecialFunctions::spherical_j_jdj_array(nmax, kr[ig], jkr[ig])
    }

    (void)wavenumber2;
    (void)wavenumber2_mr;
    (void)old_qext;
    (void)old_qsca;
    ++nmax;
  }

  return 0;
}
