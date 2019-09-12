/**
 * @file testSpecialFunctions.cpp
 *
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "SpecialFunctions.hpp"
#include <fstream>

/*! @brief Maximum order to test for Bessel and Wigner D functions,
 *  @f$n_{max}@f$. */
#define TESTSPECIALFUNCTIONS_NMAX 80u

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

  return 0;
}
