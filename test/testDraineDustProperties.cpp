/**
 * @file testDraineDustProperties.cpp
 *
 * @brief Unit test for the DraineDustProperties class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "DraineDustProperties.hpp"

#include <fstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the DraineDustProperties class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  DraineDustProperties dust_properties;
  std::ofstream dfile("test_draine_dustproperties.txt");
  for (uint_fast32_t i = 0; i < 100; ++i) {
    const float_type wavelength = 100. + 2. * (i + 0.5);
    const std::complex<float_type> m_r = dust_properties.get_refractive_index(
        wavelength * 1.e-6, 4.e-8, DUSTGRAINTYPE_CARBON);
    dfile << wavelength << "\t" << m_r.real() << "\t" << m_r.imag() << "\n";
  }

  return 0;
}
