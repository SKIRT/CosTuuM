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

  DraineHensleyShapeDistribution distribution(100);
  const float_type dmin = distribution.get_minimum_axis_ratio();
  const float_type dint = distribution.get_maximum_axis_ratio() - dmin;
  std::ofstream ofile("test_draine_hensley_shape_distribution.txt");
  ofile << "# d\tP(d)\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    const float_type d = dmin + 0.01 * (i + 0.5) * dint;
    ofile << d << "\t" << distribution(d) << "\n";
  }

  return 0;
}
