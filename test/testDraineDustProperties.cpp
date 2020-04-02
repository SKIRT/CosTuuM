/**
 * @file testDraineDustProperties.cpp
 *
 * @brief Unit test for the DraineDustProperties class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
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

  // we can test the carbon refractive index code by comparing the values
  // computed internally (based on the 0.1 micron particle table) with the
  // values in the 0.01 micron particle table
  // get reference data
  LinearInterpolatedTable<5, float_type> reference_table;
  reference_table.from_ascii_file(DRAINEDUSTPROPERTIESDATALOCATION
                                  "callindex.out_CpeD03_0.01");
  reference_table.add_column<3>(1.);

  DraineDustProperties dust_properties;
  std::ofstream dfile("test_draine_dustproperties.txt");
  dfile << "# wavelength "
           "(um)\tmr.real\tmr.imag\tmr_ref.real\tmr_ref.imag\treldiff."
           "real\treldiff.imag\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    const float_type wavelength = 10. + 2. * (i + 0.5);
    const std::complex<float_type> mr = dust_properties.get_refractive_index(
        wavelength * 1.e-6, 1.e-8, DUSTGRAINTYPE_CARBON_PERPENDICULAR);
    const TableRow<5, float_type> &refrow =
        reference_table.get_row<0>(wavelength);
    const std::complex<float_type> mr_ref(refrow[3], refrow[4]);
    const float_type reldiff_real =
        fabs(mr.real() - mr_ref.real()) / fabs(mr_ref.real());
    const float_type reldiff_imag =
        fabs(mr.imag() - mr_ref.imag()) / fabs(mr_ref.imag());
    dfile << wavelength << "\t" << mr.real() << "\t" << mr.imag() << "\t"
          << mr_ref.real() << "\t" << mr_ref.imag() << "\t" << reldiff_real
          << "\t" << reldiff_imag << "\n";
    assert_values_equal_rel(static_cast<double>(mr.real()),
                            static_cast<double>(mr_ref.real()), 0.002);
    assert_values_equal_rel(static_cast<double>(mr.imag()),
                            static_cast<double>(mr_ref.imag()), 0.002);
  }

  return 0;
}
