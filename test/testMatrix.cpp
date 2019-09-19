/**
 * @file testMatrix.cpp
 *
 * @brief Unit test for the Matrix class and its member functions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "Matrix.hpp"
#include <iostream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the Matrix class and its member functions.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Matrix<std::complex<float_type>> A(3, 3);
  A(0, 0) = 1.;
  A(0, 0).imag(1.);
  A(0, 1) = 2.;
  A(0, 1).imag(1.);
  A(0, 2) = 3.;
  A(0, 2).imag(1.);
  A(1, 0) = 4.;
  A(1, 0).imag(2.);
  A(1, 1) = 5.;
  A(1, 1).imag(1.);
  A(1, 2) = 6.;
  A(1, 2).imag(1.);
  A(2, 0) = 8.;
  A(2, 0).imag(3.);
  A(2, 1) = 9.;
  A(2, 1).imag(1.);
  A(2, 2).imag(-2.);

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  A.plu_inverse();

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
