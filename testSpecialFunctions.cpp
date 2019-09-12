/**
 * @file testSpecialFunctions.cpp
 *
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "SpecialFunctions.hpp"
#include <fstream>

/*! @brief Maximum order of Bessel function to test. The unit test only outputs
 *  this order and the first order Bessel functions and derivatives, since this
 *  covers all paths through the code. */
#define TESTSPECIALFUNCTIONS_BESSEL_NMAX 80

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
 * between the two.
 *
 * We then test the Wigner D function...
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

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
    std::complex<double> j[TESTSPECIALFUNCTIONS_BESSEL_NMAX],
        dj[TESTSPECIALFUNCTIONS_BESSEL_NMAX];
    // call the spherical Bessel function of the first kind for complex input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_BESSEL_NMAX, z,
                                            j, dj);
    // write the 1st and nth order function to the complex output file
    cfile << j[0].real() << "\t" << j[0].imag() << "\t" << dj[0].real() << "\t"
          << dj[0].imag() << "\n";
    cfile << j[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1].real() << "\t"
          << j[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1].imag() << "\t"
          << dj[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1].real() << "\t"
          << dj[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1].imag() << "\n";
    // create output arrays for the real version tests
    double jr[TESTSPECIALFUNCTIONS_BESSEL_NMAX],
        djr[TESTSPECIALFUNCTIONS_BESSEL_NMAX],
        yr[TESTSPECIALFUNCTIONS_BESSEL_NMAX],
        dyr[TESTSPECIALFUNCTIONS_BESSEL_NMAX];
    // call the spherical Bessel function of the first kind for real input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_BESSEL_NMAX, x,
                                            jr, djr);
    // call the spherical Bessel function of the second kind for real input
    // and output
    SpecialFunctions::spherical_y_ydy_array(TESTSPECIALFUNCTIONS_BESSEL_NMAX, x,
                                            yr, dyr);
    // write the 1st and nth order function to the real output file
    rfile << jr[0] << "\t" << djr[0] << "\t" << yr[0] << "\t" << dyr[0] << "\n";
    rfile << jr[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1] << "\t"
          << djr[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1] << "\t"
          << yr[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1] << "\t"
          << dyr[TESTSPECIALFUNCTIONS_BESSEL_NMAX - 1] << "\n";
  }

  return 0;
}
