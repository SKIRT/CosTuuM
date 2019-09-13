/**
 * @file testSpecialFunctions.cpp
 *
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "SpecialFunctions.hpp"
#include <fstream>

/*! @brief Maximum order to test for Bessel and Wigner D functions,
 *  @f$n_{max}@f$. */
#define TESTSPECIALFUNCTIONS_NMAX 80u

/*! @brief Order of Gauss-Legendre quadrature to test. */
#define TESTSPECIALFUNCTIONS_NGAUSS 100u

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
    const double x = 0.05 + 0.1 * i;
    // we create a complex argument to test the complex version
    const std::complex<double> z(x, x);
    // create output arrays for the complex version test
    std::complex<double> j[TESTSPECIALFUNCTIONS_NMAX],
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
    double jr[TESTSPECIALFUNCTIONS_NMAX], djr[TESTSPECIALFUNCTIONS_NMAX],
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
    const double cosx = -0.999 + 0.002 * i;
    // create output arrays
    double d[TESTSPECIALFUNCTIONS_NMAX], dd[TESTSPECIALFUNCTIONS_NMAX];
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

  /// Equal volume sphere to equal surface area sphere radius ratio

  assert_condition(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          1.) == 1.);
  assert_values_equal_rel(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          0.5),
      0.9637112829756893, 1.e-10);
  assert_values_equal_rel(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(2.),
      0.955443263377723, 1.e-10);

  /// Gauss-Legendre quadrature

  std::vector<double> x(TESTSPECIALFUNCTIONS_NGAUSS),
      w(TESTSPECIALFUNCTIONS_NGAUSS);
  // first test against the values on Wikipedia
  {
    // n = 1
    SpecialFunctions::get_gauss_legendre_points_and_weigths(1, x, w);
    assert_condition(x[0] == 0.);
    assert_condition(w[0] == 2.);
  }
  {
    // n = 2
    SpecialFunctions::get_gauss_legendre_points_and_weigths(2, x, w);
    assert_values_equal_rel(x[0], -1. / std::sqrt(3.), 1.e-10);
    assert_values_equal_rel(w[0], 1., 1.e-10);
    assert_values_equal_rel(x[1], 1. / std::sqrt(3.), 1.e-10);
    assert_values_equal_rel(w[1], 1., 1.e-10);
  }
  {
    // n = 3
    SpecialFunctions::get_gauss_legendre_points_and_weigths(3, x, w);
    assert_values_equal_rel(x[0], -std::sqrt(3. / 5.), 1.e-10);
    assert_values_equal_rel(w[0], 5. / 9., 1.e-10);
    assert_condition(x[1] == 0.);
    assert_values_equal_rel(w[1], 8. / 9., 1.e-10);
    assert_values_equal_rel(x[2], std::sqrt(3. / 5.), 1.e-10);
    assert_values_equal_rel(w[2], 5. / 9., 1.e-10);
  }
  {
    // n = 4
    SpecialFunctions::get_gauss_legendre_points_and_weigths(4, x, w);
    assert_values_equal_rel(
        x[0], -std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)), 1.e-10);
    assert_values_equal_rel(w[0], (18 - std::sqrt(30.)) / 36., 1.e-10);
    assert_values_equal_rel(
        x[1], -std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)), 1.e-10);
    assert_values_equal_rel(w[1], (18. + std::sqrt(30.)) / 36., 1.e-10);
    assert_values_equal_rel(
        x[2], std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)), 1.e-10);
    assert_values_equal_rel(w[2], (18. + std::sqrt(30.)) / 36., 1.e-10);
    assert_values_equal_rel(
        x[3], std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)), 1.e-10);
    assert_values_equal_rel(w[3], (18 - std::sqrt(30.)) / 36., 1.e-10);
  }
  {
    // n = 5
    SpecialFunctions::get_gauss_legendre_points_and_weigths(5, x, w);
    assert_values_equal_rel(
        x[0], -1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)), 1.e-10);
    assert_values_equal_rel(w[0], (322. - 13. * std::sqrt(70.)) / 900., 1.e-10);
    assert_values_equal_rel(
        x[1], -1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)), 1.e-10);
    assert_values_equal_rel(w[1], (322. + 13. * std::sqrt(70.)) / 900., 1.e-10);
    assert_condition(x[2] == 0.);
    assert_values_equal_rel(w[2], 128. / 225., 1.e-10);
    assert_values_equal_rel(
        x[3], 1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)), 1.e-10);
    assert_values_equal_rel(w[3], (322. + 13. * std::sqrt(70.)) / 900., 1.e-10);
    assert_values_equal_rel(
        x[4], 1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)), 1.e-10);
    assert_values_equal_rel(w[4], (322. - 13. * std::sqrt(70.)) / 900., 1.e-10);
  }

  // now print additional values for a much higher order to analyse with the
  // script
  SpecialFunctions::get_gauss_legendre_points_and_weigths(
      TESTSPECIALFUNCTIONS_NGAUSS, x, w);
  std::ofstream gfile("test_gauss_legendre_quadrature.txt");
  for (uint_fast32_t i = 0; i < TESTSPECIALFUNCTIONS_NGAUSS; ++i) {
    gfile << x[i] << "\t" << w[i] << "\n";
  }

  return 0;
}
