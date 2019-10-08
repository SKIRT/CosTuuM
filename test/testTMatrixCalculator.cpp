/**
 * @file testTMatrixCalculator.cpp
 *
 * @brief Unit test for the TMatrixCalculator class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "OrientationDistribution.hpp"
#include "TMatrixCalculator.hpp"
#include "UnitConverter.hpp"

#include <fstream>
#include <sstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the TMatrixCalculator class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  std::ifstream ifile("test_tmatrixcalculator.txt");
  std::string line;
  // skip the first comment line
  getline(ifile, line);
  // now read the other lines
  uint_fast32_t counter = 0;
  while (getline(ifile, line)) {
    ++counter;

    ctm_warning("Line %" PRIuFAST32, counter);

    std::istringstream linestream(line);
    float_type axi, rat, lam, mrr, mri, eps, ddelt, alpha, beta, thet0, thet,
        phi0, phi, refqsca, refqext, refwalb, refZ[4][4];
    uint_fast32_t ndgs;
    linestream >> axi >> rat >> lam >> mrr >> mri >> eps >> ddelt >> ndgs >>
        alpha >> beta >> thet0 >> thet >> phi0 >> phi >> refqsca >> refqext >>
        refwalb >> refZ[0][0] >> refZ[0][1] >> refZ[0][2] >> refZ[0][3] >>
        refZ[1][0] >> refZ[1][1] >> refZ[1][2] >> refZ[1][3] >> refZ[2][0] >>
        refZ[2][1] >> refZ[2][2] >> refZ[2][3] >> refZ[3][0] >> refZ[3][1] >>
        refZ[3][2] >> refZ[3][3];

    TMatrix *Tmat = TMatrixCalculator::calculate_TMatrix(
        rat, eps, axi, lam, 200, ddelt * 0.1, ndgs,
        std::complex<float_type>(mrr, mri), 500);

    const float_type qext = Tmat->get_extinction_coefficient();
    const float_type qsca = Tmat->get_scattering_coefficient();
    const float_type walb = -qsca / qext;

    alpha = UnitConverter::to_SI<QUANTITY_ANGLE>(double(alpha), "degrees");
    beta = UnitConverter::to_SI<QUANTITY_ANGLE>(double(beta), "degrees");
    thet0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet0), "degrees");
    thet = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet), "degrees");
    phi0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi0), "degrees");
    phi = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi), "degrees");

    assert_values_equal_rel(double(qext), double(refqext), 1.e-5);
    assert_values_equal_rel(double(qsca), double(refqsca), 1.e-5);
    assert_values_equal_rel(double(walb), double(refwalb), 1.e-5);

    Matrix<float_type> Z =
        Tmat->get_scattering_matrix(alpha, beta, thet0, phi0, thet, phi);

    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        assert_values_equal_rel(double(Z(i, j)), double(refZ[i][j]), 2.e-2);
      }
    }

    if (counter == 1) {
      OrientationDistribution orientation_distribution(2 * Tmat->get_nmax());
      TMatrix *T_ensemble = TMatrixCalculator::apply_orientation_distribution(
          *Tmat, orientation_distribution);
      delete T_ensemble;
    }

    delete Tmat;
  }

  return 0;
}
