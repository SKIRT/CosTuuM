/**
 * @file CTM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Matrix.hpp"
#include "SpecialFunctions.hpp"
#include "TMatrix.hpp"

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
  const uint_fast32_t ndgs = 2;

  /// hardcoded program parameters
  // maximum number of iterations during the first loop
  const uint_fast32_t maximum_order = 200;
  // maximum number of iterations during the second loop
  const uint_fast32_t maximum_ngauss = 500;

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

    TMatrix T(wavelength, mr, R_V, axis_ratio, nmax, ngauss);

    double qsca = 0.;
    double qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const double dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    dsca = std::abs((old_qsca - qsca) / qsca);
    dext = std::abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    ctm_warning("dsca: %g, dext: %g", dsca, dext);

    ++nmax;
  }

  if (nmax == maximum_order) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in last iteration
    --nmax;
    ctm_warning("Converged for nmax = %" PRIuFAST32, nmax);
  }

  uint_fast32_t ngauss = nmax * ndgs + 1;
  dext = tolerance + 1.;
  dsca = tolerance + 1.;
  while (ngauss < maximum_ngauss && (dext > tolerance || dsca > tolerance)) {

    TMatrix T(wavelength, mr, R_V, axis_ratio, nmax, ngauss);

    double qsca = 0.;
    double qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const double dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    dsca = std::abs((old_qsca - qsca) / qsca);
    dext = std::abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    ctm_warning("dsca: %g, dext: %g", dsca, dext);

    ++ngauss;
  }

  if (ngauss == maximum_ngauss) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in final iteration
    --ngauss;
    ctm_warning("Converged for ngauss = %" PRIuFAST32, ngauss);
  }

  return 0;
}
