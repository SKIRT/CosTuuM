/**
 * @file CTM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include <cinttypes>
#include <cmath>
#include <complex>

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
  double rat = 0.1;
  // wavelength of incoming radiation (in same units as the particle size)
  const double lam = 2. * M_PI;
  // refractive index
  const std::complex<double> mr(1.5, 0.02);
  // ratio of horizontal and vertical axis of the spheroidal particle
  const double eps = 0.5;
  // tolerance for the calculation
  const double ddelt = 1.e-4;
  // number of Gauss-Legendre points to use as a multiplier of the maximum
  // order of spherical harmonics used to decompose the electric field, during
  // the first loop of the algorithm
  const uint_fast32_t ndgs = 1;

  /// hardcoded program parameters
  // maximum number of iterations during the first loop
  const uint_fast32_t npn1 = 200;
  // maximum number of iterations during the second loop
  const uint_fast32_t npng1 = 500;

  // make sure 'rat' contains the right ratio if it is not 1
  if (std::abs(rat - 1.) > 1.e-8) {
    //    rat =
    //    SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(eps);
  }

  return 0;
}
