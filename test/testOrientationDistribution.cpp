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

  std::ofstream ofile("test_orientationdistribution.txt");
  ofile << "# n\tp_n\n";
  for (uint_fast32_t i = 0; i < 101; ++i) {
    ofile << i << "\t" << od.get_coefficient(i) << "\n";
  }

  return 0;
}
