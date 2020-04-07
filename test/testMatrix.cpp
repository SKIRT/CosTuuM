/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

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

  Matrix<std::complex<float_type>> A(4, 4);
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

  // add some additional elements to the part of the matrix we don't care about
  A(3, 0) = 42.;
  A(3, 1) = 42.;
  A(3, 2) = 48.;
  A(3, 3) = 44.;
  A(0, 3) = 41.;
  A(1, 3) = 48.;
  A(2, 3) = 43.;

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  A.plu_inverse(3);

  A.binary_dump("test_matrix.dump");

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
