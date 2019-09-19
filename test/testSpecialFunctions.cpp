/**
 * @file testSpecialFunctions.cpp
 *
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "Configuration.hpp"
#include "SpecialFunctions.hpp"
#include <fstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/*! @brief Maximum order to test for Bessel and Wigner D functions,
 *  @f$n_{max}@f$. */
const uint_fast32_t TESTSPECIALFUNCTIONS_NMAX = 80;

/*! @brief Order of Gauss-Legendre quadrature to test. */
const uint_fast32_t TESTSPECIALFUNCTIONS_NGAUSS = 200;

/**
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * We call the functions SpecialFunctions::spherical_j_jdj_array() and
 * SpecialFunctions::spherical_y_ydy_array() for a range of input values (real
 * and complex for the first kind), and write them to two text files, called
 * test_bessel_complex.txt and test_bessel_real.txt. The accompanying script
 * plot_test_special_functions.py reads these files and computes the same
 * functions for the same input values using the Bessel functions that are part
 * of scipy.special. It then plots both versions and the relative difference
 * between the two. Note that we need to test both the first order and
 * @f$n_max@f$th order derivatives, as they use different bits of code.
 *
 * We then test the Wigner D function by calling
 * SpecialFunctions::wigner_dn_0m() on a range of input values, and write the
 * functions and their derivatives to a text file called test_wigner_d.txt.
 * The accompanying script reads this file and computes the same functions for
 * the same input values using the associated Legendre polynomials in
 * scipy.special and the relation between Wigner D functions and associated
 * Legendre polynomials. It then plots both versions and the relative
 * difference between the two. We test both @f$m = 0@f$ and @f$m = n_{max}@f$
 * to cover all possible paths through the code.
 *
 * We then test
 * SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio() by
 * calling it with 3 different axis ratios (1, <1, >1) and comparing with
 * results obtained with Python.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// Bessel functions

  // open the output text files
  std::ofstream cfile("test_bessel_complex.txt");
  std::ofstream rfile("test_bessel_real.txt");
  // loop over a spatial range
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(0.05, 100., 0.1) in the
    // Python script
    const float_type x = 0.05 + 0.1 * i;
    // we create a complex argument to test the complex version
    const std::complex<float_type> z(x, x);
    // create output arrays for the complex version test
    std::complex<float_type> j[TESTSPECIALFUNCTIONS_NMAX],
        dj[TESTSPECIALFUNCTIONS_NMAX];
    // call the spherical Bessel function of the first kind for complex input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_NMAX, z, j,
                                            dj);
    // write the 1st and nth order function to the complex output file
    cfile << j[0].real() << "\t" << j[0].imag() << "\t" << dj[0].real() << "\t"
          << dj[0].imag() << "\n";
    cfile << j[TESTSPECIALFUNCTIONS_NMAX - 1].real() << "\t"
          << j[TESTSPECIALFUNCTIONS_NMAX - 1].imag() << "\t"
          << dj[TESTSPECIALFUNCTIONS_NMAX - 1].real() << "\t"
          << dj[TESTSPECIALFUNCTIONS_NMAX - 1].imag() << "\n";
    // create output arrays for the real version tests
    float_type jr[TESTSPECIALFUNCTIONS_NMAX], djr[TESTSPECIALFUNCTIONS_NMAX],
        yr[TESTSPECIALFUNCTIONS_NMAX], dyr[TESTSPECIALFUNCTIONS_NMAX];
    // call the spherical Bessel function of the first kind for real input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_NMAX, x, jr,
                                            djr);
    // call the spherical Bessel function of the second kind for real input
    // and output
    SpecialFunctions::spherical_y_ydy_array(TESTSPECIALFUNCTIONS_NMAX, x, yr,
                                            dyr);
    // write the 1st and nth order function to the real output file
    rfile << jr[0] << "\t" << djr[0] << "\t" << yr[0] << "\t" << dyr[0] << "\n";
    rfile << jr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << djr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << yr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dyr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  /// Wigner D function

  // open the output text file
  std::ofstream dfile("test_wigner_d.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(-0.999, 1., 0.002) expression
    // in the Python script
    const float_type cosx = -0.999 + 0.002 * i;
    // create output arrays
    float_type d[TESTSPECIALFUNCTIONS_NMAX], dd[TESTSPECIALFUNCTIONS_NMAX];
    // call the function with m=0
    SpecialFunctions::wigner_dn_0m(cosx, TESTSPECIALFUNCTIONS_NMAX, 0, d, dd);
    // write an output line
    dfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
    // call the function with m=nmax
    SpecialFunctions::wigner_dn_0m(cosx, TESTSPECIALFUNCTIONS_NMAX,
                                   TESTSPECIALFUNCTIONS_NMAX, d, dd);
    // write an output line
    dfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  // open the output text file
  std::ofstream dsfile("test_wigner_d_sinx.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(-0.999, 1., 0.002) expression
    // in the Python script
    const float_type cosx = -0.999 + 0.002 * i;
    // create output arrays
    float_type d[TESTSPECIALFUNCTIONS_NMAX], dd[TESTSPECIALFUNCTIONS_NMAX];
    // call the function with m=0
    SpecialFunctions::wigner_dn_0m_sinx(cosx, TESTSPECIALFUNCTIONS_NMAX, 0, d,
                                        dd);
    // write an output line
    dsfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
           << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
    // call the function with m=nmax
    SpecialFunctions::wigner_dn_0m_sinx(cosx, TESTSPECIALFUNCTIONS_NMAX,
                                        TESTSPECIALFUNCTIONS_NMAX, d, dd);
    // write an output line
    dsfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
           << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  /// Equal volume sphere to equal surface area sphere radius ratio

  assert_condition(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(1.)) == float_type(1.));
  const double prolate = double(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(0.5)));
  assert_values_equal_rel(prolate, 0.9637112829756893, 1.e-10);
  const double oblate = double(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(2.)));
  assert_values_equal_rel(oblate, 0.955443263377723, 1.e-10);

  /// Gauss-Legendre quadrature

  const double tolerance = 1.e-7;
  std::vector<float_type> x(TESTSPECIALFUNCTIONS_NGAUSS),
      w(TESTSPECIALFUNCTIONS_NGAUSS);
  // first test against the values on Wikipedia
  {
    // n = 1
    SpecialFunctions::get_gauss_legendre_points_and_weights(1, x, w);
    ctm_warning("x: [%g], w: [%g]", double(x[0]), double(w[0]));
    assert_condition(x[0] == 0.);
    assert_condition(w[0] == 2.);
  }
  {
    // n = 2
    SpecialFunctions::get_gauss_legendre_points_and_weights(2, x, w);
    ctm_warning("x: [%g, %g], w: [%g, %g]", double(x[0]), double(x[1]),
                double(w[0]), double(w[1]));
    assert_values_equal_rel(double(x[0]), -1. / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[0]), 1., tolerance);
    assert_values_equal_rel(double(x[1]), 1. / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[1]), 1., tolerance);
  }
  {
    // n = 3
    SpecialFunctions::get_gauss_legendre_points_and_weights(3, x, w);
    ctm_warning("x: [%g, %g, %g], w: [%g, %g, %g]", double(x[0]), double(x[1]),
                double(x[2]), double(w[0]), double(w[1]), double(w[2]));
    assert_values_equal_rel(double(x[0]), -std::sqrt(3. / 5.), tolerance);
    assert_values_equal_rel(double(w[0]), 5. / 9., tolerance);
    assert_values_equal_tol(double(x[1]), 0., tolerance);
    assert_values_equal_rel(double(w[1]), 8. / 9., tolerance);
    assert_values_equal_rel(double(x[2]), std::sqrt(3. / 5.), tolerance);
    assert_values_equal_rel(double(w[2]), 5. / 9., tolerance);
  }
  {
    // n = 4
    SpecialFunctions::get_gauss_legendre_points_and_weights(4, x, w);
    ctm_warning("x: [%g, %g, %g, %g], w: [%g, %g, %g, %g]", double(x[0]),
                double(x[1]), double(x[2]), double(x[3]), double(w[0]),
                double(w[1]), double(w[2]), double(w[3]));
    assert_values_equal_rel(double(x[0]),
                            -std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[0]), (18 - std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[1]),
                            -std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[1]), (18. + std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[2]),
                            std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[2]), (18. + std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[3]),
                            std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[3]), (18 - std::sqrt(30.)) / 36.,
                            tolerance);
  }
  {
    // n = 5
    SpecialFunctions::get_gauss_legendre_points_and_weights(5, x, w);
    ctm_warning("x: [%g, %g, %g, %g, %g], w: [%g, %g, %g, %g, %g]",
                double(x[0]), double(x[1]), double(x[2]), double(x[3]),
                double(x[4]), double(w[0]), double(w[1]), double(w[2]),
                double(w[3]), double(w[4]));
    assert_values_equal_rel(double(x[0]),
                            -1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[0]), (322. - 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_rel(double(x[1]),
                            -1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[1]), (322. + 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_tol(double(x[2]), 0., tolerance);
    assert_values_equal_rel(double(w[2]), 128. / 225., tolerance);
    assert_values_equal_rel(double(x[3]),
                            1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[3]), (322. + 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_rel(double(x[4]),
                            1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[4]), (322. - 13. * std::sqrt(70.)) / 900.,
                            tolerance);
  }

  // now print additional values for a much higher order to analyse with the
  // script
  SpecialFunctions::get_gauss_legendre_points_and_weights(
      TESTSPECIALFUNCTIONS_NGAUSS, x, w);
  std::ofstream gfile("test_gauss_legendre_quadrature.txt");
  for (uint_fast32_t i = 0; i < TESTSPECIALFUNCTIONS_NGAUSS; ++i) {
    gfile << x[i] << "\t" << w[i] << "\n";
  }

  return 0;
}
