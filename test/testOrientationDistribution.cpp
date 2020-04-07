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
 * @file testOrientationDistribution.cpp
 *
 * @brief Unit test for the OrientationDistribution class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "OrientationDistribution.hpp"

#include <fstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the OrientationDistribution class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  OrientationDistribution od(100);
  od.initialise();

  std::ofstream ofile("test_orientationdistribution.txt");
  ofile << "# n\tp_n\n";
  for (uint_fast32_t i = 0; i < 101; ++i) {
    ofile << i << "\t" << od.get_coefficient(i) << "\n";
  }

  std::ofstream rfile("test_orientationdistribution_ref.txt");
  rfile << "# beta\tp(beta)\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const float_type beta = 0.001 * (i + 0.5) * M_PI;
    rfile << beta << "\t" << od(beta) << "\n";
  }

  return 0;
}
