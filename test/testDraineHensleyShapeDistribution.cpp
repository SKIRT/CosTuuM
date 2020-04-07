/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDraineHensleyShapeDistribution.cpp
 *
 * @brief Unit test for the DraineHensleyShapeDistribution class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "DraineHensleyShapeDistribution.hpp"

#include <fstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the DraineHensleyShapeDistribution class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // test the default cutoff scheme
  {
    DraineHensleyShapeDistribution distribution(100);
    const float_type dmin = distribution.get_minimum_axis_ratio();
    const float_type dint = distribution.get_maximum_axis_ratio() - dmin;
    std::ofstream ofile("test_draine_hensley_shape_distribution.txt");
    ofile << "# d\tP(d)\n";
    for (uint_fast32_t i = 0; i < 100; ++i) {
      const float_type d = dmin + 0.01 * (i + 0.5) * dint;
      ofile << d << "\t" << distribution(d) << "\n";
    }
  }

  // test the advanced fraction scheme
  {
    DraineHensleyShapeDistribution distribution(100, 0., 0.9);
    const float_type dmin = distribution.get_minimum_axis_ratio();
    const float_type dint = distribution.get_maximum_axis_ratio() - dmin;
    std::ofstream ofile("test_draine_hensley_shape_distribution_fraction.txt");
    ofile << "# d\tP(d)\n";
    for (uint_fast32_t i = 0; i < 100; ++i) {
      const float_type d = dmin + 0.01 * (i + 0.5) * dint;
      ofile << d << "\t" << distribution(d) << "\n";
    }
  }

  return 0;
}
