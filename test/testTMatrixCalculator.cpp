/**
 * @file testTMatrixCalculator.cpp
 *
 * @brief Unit test for the TMatrixCalculator class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "TMatrixCalculator.hpp"

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
  while (getline(ifile, line)) {
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
        rat, eps, axi, lam, 200, ddelt, ndgs,
        std::complex<float_type>(mrr, mri), 500);

    const float_type qext = Tmat->get_extinction_coefficient();
    const float_type qsca = Tmat->get_scattering_coefficient();
    const float_type walb = -qsca / qext;

    const float_type deg_to_rad = M_PI / 180.;
    alpha *= deg_to_rad;
    beta *= deg_to_rad;
    thet0 *= deg_to_rad;
    thet *= deg_to_rad;
    phi0 *= deg_to_rad;
    phi *= deg_to_rad;

    assert_values_equal_rel(double(qext), double(refqext), 1.e-3);
    assert_values_equal_rel(double(qsca), double(refqsca), 1.e-3);
    assert_values_equal_rel(double(walb), double(refwalb), 1.e-3);

    Matrix<float_type> Z =
        Tmat->get_scattering_matrix(alpha, beta, thet0, phi0, thet, phi);

    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        assert_values_equal_rel(double(Z(i, j)), double(refZ[i][j]), 0.1);
      }
    }

    delete Tmat;
  }

  return 0;
}
