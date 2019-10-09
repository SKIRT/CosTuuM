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

  /// output forward scattering coefficients for a single and ensemble T-matrix
  {
    TMatrix *Tmatrix_single = TMatrixCalculator::calculate_TMatrix(
        0.1, 0.5, 10., 2. * M_PI, 200, 1.e-4, 2,
        std::complex<float_type>(1.5, 0.02), 500);
    OrientationDistribution orientation_distribution(
        2 * Tmatrix_single->get_nmax());
    TMatrix *Tmatrix_ensemble =
        TMatrixCalculator::apply_orientation_distribution(
            *Tmatrix_single, orientation_distribution);

    // check backscattering matrix property
    {
      Matrix<float_type> Zsingle = Tmatrix_single->get_scattering_matrix(
          0., 0., 0.5 * M_PI, 0., 0.5 * M_PI, M_PI);
      const float_type Zsingle_sum =
          Zsingle(0, 0) - Zsingle(1, 1) + Zsingle(2, 2) - Zsingle(3, 3);
      assert_values_equal_tol(double(Zsingle_sum), 0., 1.e-10);
    }
    {
      Matrix<float_type> Zensemble = Tmatrix_ensemble->get_scattering_matrix(
          0., 0., 0.5 * M_PI, 0., 0.5 * M_PI, M_PI);
      const float_type Zsingle_sum =
          Zensemble(0, 0) - Zensemble(1, 1) + Zensemble(2, 2) - Zensemble(3, 3);
      assert_values_equal_tol(double(Zsingle_sum), 0., 1.e-10);
    }

    // check that extinction matrix satisfies theoretical criterion
    {
      Matrix<float_type> K =
          Tmatrix_ensemble->get_extinction_matrix(0., 0., 0.5 * M_PI, 0.);
      assert_values_equal_tol(double(K(0, 2)), 0., 1.e-10);
      assert_values_equal_tol(double(K(0, 3)), 0., 1.e-10);
      assert_values_equal_tol(double(K(1, 2)), 0., 1.e-10);
      assert_values_equal_tol(double(K(1, 3)), 0., 1.e-10);
      assert_values_equal_tol(double(K(2, 0)), 0., 1.e-10);
      assert_values_equal_tol(double(K(3, 0)), 0., 1.e-10);
      assert_values_equal_tol(double(K(2, 1)), 0., 1.e-10);
      assert_values_equal_tol(double(K(3, 1)), 0., 1.e-10);
    }

    std::vector<float_type> xphi(100), wphi(100);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        100, 0., 2. * M_PI, xphi, wphi);
    std::vector<float_type> xtheta(100), wtheta(100);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        100, 0., M_PI, xtheta, wtheta);
    float_type Zsinglequad = 0.;
    float_type Zensemblequad = 0.;
    std::ofstream ofile("test_tmatrixcalculator_result.txt");
    ofile << "# theta\tphi\tZ00\n";
    for (uint_fast32_t i = 0; i < 100; ++i) {
      const float_type theta_out = xtheta[i];
      const float_type sintheta = sin(theta_out);
      for (uint_fast32_t j = 0; j < 100; ++j) {
        const float_type phi_out = xphi[j];
        Matrix<float_type> Zsingle = Tmatrix_single->get_scattering_matrix(
            0., 0., 0.5 * M_PI, 0., theta_out, phi_out);
        Matrix<float_type> Zensemble = Tmatrix_ensemble->get_scattering_matrix(
            0., 0., 0.5 * M_PI, 0., theta_out, phi_out);
        ofile << theta_out << "\t" << phi_out << "\t" << Zsingle(0, 0) << "\t"
              << Zensemble(0, 0) << "\n";
        Zsinglequad += wphi[j] * wtheta[i] * sintheta * Zsingle(0, 0);
        Zensemblequad += wphi[j] * wtheta[i] * sintheta * Zensemble(0, 0);
      }
    }
    ctm_warning("Quad: %g %g", double(Zsinglequad), double(Zensemblequad));

    delete Tmatrix_single;
    delete Tmatrix_ensemble;
  }

  return 0;
}
