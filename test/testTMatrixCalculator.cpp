/**
 * @file testTMatrixCalculator.cpp
 *
 * @brief Unit test for the TMatrixCalculator class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "Configuration.hpp"
#include "OrientationDistribution.hpp"
#include "ShapeDistribution.hpp"
#include "TMatrixCalculator.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"

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

  // control which tests are performed
  // ideally, this is all of them, but sometimes (temporarily) disabling some
  // can be useful
  const bool do_benchmark_test = true;
  const bool do_ensemble_test = true;
  const bool do_astro_test = true;
  const bool do_shape_test = true;
  const bool do_fixed_orientation_test = true;

  /// benchmark test: check T matrix results against some results obtained
  /// with Mishchenko's original T matrix code for the same input values
  if (do_benchmark_test) {
    // file containing the benchmark test data
    std::ifstream ifile("test_tmatrixcalculator.txt");
    std::string line;
    // skip the first comment line
    getline(ifile, line);
    // now read the other lines
    uint_fast32_t counter = 0;
    while (getline(ifile, line)) {
      ++counter;

      ctm_warning("Line %" PRIuFAST32, counter);

      // read the input parameters and reference output
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

      // convert the input parameters to SI units
      axi = UnitConverter::to_SI<QUANTITY_LENGTH>(double(axi), "micron");
      lam = UnitConverter::to_SI<QUANTITY_LENGTH>(double(lam), "micron");
      // convert the reference results to SI units
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          refZ[i][j] = UnitConverter::to_SI<QUANTITY_SURFACE_AREA>(
              double(refZ[i][j]), "micron^2");
        }
      }

      // construct the T matrix for the same input parameters
      const TMatrix *Tmat = TMatrixCalculator::calculate_TMatrix(
          rat, eps, axi, lam, 200, ddelt * 0.1, ndgs,
          std::complex<float_type>(mrr, mri), 500);

      // check that the extinction and scattering coefficients match the
      // expected values
      const float_type qext = Tmat->get_extinction_coefficient();
      const float_type qsca = Tmat->get_scattering_coefficient();
      const float_type walb = -qsca / qext;

      assert_values_equal_rel(double(qext), double(refqext), 1.e-5);
      assert_values_equal_rel(double(qsca), double(refqsca), 1.e-5);
      assert_values_equal_rel(double(walb), double(refwalb), 1.e-5);

      // convert the input parameter angles to radians
      alpha = UnitConverter::to_SI<QUANTITY_ANGLE>(double(alpha), "degrees");
      beta = UnitConverter::to_SI<QUANTITY_ANGLE>(double(beta), "degrees");
      thet0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet0), "degrees");
      thet = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet), "degrees");
      phi0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi0), "degrees");
      phi = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi), "degrees");

      // compute the scattering matrix
      const Matrix<float_type> Z =
          Tmat->get_scattering_matrix(alpha, beta, thet0, phi0, thet, phi);

      // compare the result with the reference
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          assert_values_equal_rel(double(Z(i, j)), double(refZ[i][j]), 2.e-2);
        }
      }

      delete Tmat;
    }
  }

  /// output forward scattering coefficients for a single and ensemble T-matrix
  if (do_ensemble_test) {
    ctm_warning("Calculating initial T matrix...");
    const TMatrix *Tmatrix_single = TMatrixCalculator::calculate_TMatrix(
        0.1, 0.5, 1.e-5, 2.e-6 * M_PI, 200, 1.e-4, 2,
        std::complex<float_type>(1.5, 0.02), 500);
    ctm_warning("Done.");

    ctm_warning("Averaging T matrix over orientations...");
    OrientationDistribution orientation_distribution(
        2 * Tmatrix_single->get_nmax());
    orientation_distribution.initialise();
    const TMatrix *Tmatrix_ensemble =
        TMatrixCalculator::apply_orientation_distribution(
            *Tmatrix_single, orientation_distribution);
    ctm_warning("Done.");

    ctm_warning("Computing ensemble extinction matrix...");
    const Matrix<float_type> Kensemble =
        Tmatrix_ensemble->get_extinction_matrix(0., 0., 0.3 * M_PI, 0.);
    ctm_warning("Done.");

    ctm_warning("Computing reference averaged extinction matrix...")
        const uint_fast32_t ngauss_beta = 100;
    const uint_fast32_t ngauss_alpha = 200;
    Matrix<float_type> Kref(4, 4);
    std::vector<float_type> beta(ngauss_beta), weights_beta(ngauss_beta);
    std::vector<float_type> alpha(ngauss_alpha), weights_alpha(ngauss_alpha);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        ngauss_beta, 0., M_PI, beta, weights_beta);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        ngauss_alpha, 0., 2. * M_PI, alpha, weights_alpha);
    for (uint_fast32_t iga = 0; iga < ngauss_alpha; ++iga) {
      for (uint_fast32_t igb = 0; igb < ngauss_beta; ++igb) {
        const Matrix<float_type> Ksingle =
            Tmatrix_single->get_extinction_matrix(alpha[iga], beta[igb],
                                                  0.3 * M_PI, 0.);
        const float_type factor = 0.5 * M_1_PI * sin(beta[igb]) *
                                  orientation_distribution(beta[igb]) *
                                  weights_beta[igb] * weights_alpha[iga];
        for (uint_fast8_t i = 0; i < 4; ++i) {
          for (uint_fast8_t j = 0; j < 4; ++j) {
            Kref(i, j) += factor * Ksingle(i, j);
          }
        }
      }
    }
    ctm_warning("Done.");

    ctm_warning("Result:");
    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        ctm_warning("K(%" PRIuFAST8 ",%" PRIuFAST8 "): %g %g", i, j,
                    double(Kensemble(i, j)), double(Kref(i, j)));
        assert_values_equal_tol(double(Kensemble(i, j)), double(Kref(i, j)),
                                1.e-6);
      }
    }

    // check backscattering matrix property
    {
      const Matrix<float_type> Zsingle = Tmatrix_single->get_scattering_matrix(
          0., 0., 0.5 * M_PI, 0., 0.5 * M_PI, M_PI);
      const float_type Zsingle_sum =
          Zsingle(0, 0) - Zsingle(1, 1) + Zsingle(2, 2) - Zsingle(3, 3);
      assert_values_equal_tol(double(Zsingle_sum), 0., 1.e-10);
    }
    {
      const Matrix<float_type> Zensemble =
          Tmatrix_ensemble->get_scattering_matrix(0., 0., 0.5 * M_PI, 0.,
                                                  0.5 * M_PI, M_PI);
      const float_type Zsingle_sum =
          Zensemble(0, 0) - Zensemble(1, 1) + Zensemble(2, 2) - Zensemble(3, 3);
      assert_values_equal_tol(double(Zsingle_sum), 0., 1.e-10);
    }

    // check that extinction matrix satisfies theoretical criterion
    {
      const Matrix<float_type> K =
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
    std::ofstream ofile("test_tmatrixcalculator_result.txt");
    ofile << "# theta\tphi\tZ00\n";
    for (uint_fast32_t i = 0; i < 100; ++i) {
      const float_type theta_out = xtheta[i];
      const float_type sintheta = sin(theta_out);
      for (uint_fast32_t j = 0; j < 100; ++j) {
        const float_type phi_out = xphi[j];
        const Matrix<float_type> Zsingle =
            Tmatrix_single->get_scattering_matrix(0., 0., 0.5 * M_PI, 0.,
                                                  theta_out, phi_out);
        const Matrix<float_type> Ksingle =
            Tmatrix_single->get_extinction_matrix(0., 0., theta_out, phi_out);
        ofile << theta_out << "\t" << phi_out;
        for (uint_fast8_t irow = 0; irow < 4; ++irow) {
          for (uint_fast8_t icol = 0; icol < 4; ++icol) {
            ofile << "\t" << Zsingle(irow, icol);
          }
        }
        for (uint_fast8_t irow = 0; irow < 4; ++irow) {
          for (uint_fast8_t icol = 0; icol < 4; ++icol) {
            ofile << "\t" << Ksingle(irow, icol);
          }
        }
        ofile << "\n";
        Zsinglequad += wphi[j] * wtheta[i] * sintheta * Zsingle(0, 0);
      }
    }
    ctm_warning("Quad: %g", double(Zsinglequad));

    delete Tmatrix_single;
    delete Tmatrix_ensemble;
  }

  /// astrophysically relevant T-matrix
  if (do_astro_test) {
    const TMatrix *Tmatrix = TMatrixCalculator::calculate_TMatrix(
        1., 0.5, 2.e-7, 1.e-4, 200, 1.e-4, 2, std::complex<float_type>(4., 0.1),
        500);
    std::vector<float_type> xphi(100), wphi(100);

    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        100, 0., 2. * M_PI, xphi, wphi);
    std::vector<float_type> xtheta(100), wtheta(100);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        100, 0., M_PI, xtheta, wtheta);
    std::ofstream ofile("test_tmatrixcalculator_ISM.txt");
    ofile << "# theta\tphi\tZ\tK\n";
    for (uint_fast32_t i = 0; i < 100; ++i) {
      const float_type theta_out = xtheta[i];
      for (uint_fast32_t j = 0; j < 100; ++j) {
        const float_type phi_out = xphi[j];
        const Matrix<float_type> Zsingle = Tmatrix->get_scattering_matrix(
            0., 0., 0.5 * M_PI, 0., theta_out, phi_out);
        const Matrix<float_type> Ksingle =
            Tmatrix->get_extinction_matrix(0., 0., theta_out, phi_out);
        ofile << theta_out << "\t" << phi_out;
        for (uint_fast8_t irow = 0; irow < 4; ++irow) {
          for (uint_fast8_t icol = 0; icol < 4; ++icol) {
            ofile << "\t" << Zsingle(irow, icol);
          }
        }
        for (uint_fast8_t irow = 0; irow < 4; ++irow) {
          for (uint_fast8_t icol = 0; icol < 4; ++icol) {
            ofile << "\t" << Ksingle(irow, icol);
          }
        }
        ofile << "\n";
      }
    }
    delete Tmatrix;
  }

  /// test the procedure to average over a shape distribution
  if (do_shape_test) {

    ShapeDistribution shape_distribution;

    const uint_fast32_t ngauss = 100;
    std::vector<float_type> axis_ratio(ngauss), weights(ngauss);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        ngauss, shape_distribution.get_minimum_axis_ratio(),
        shape_distribution.get_maximum_axis_ratio(), axis_ratio, weights);
    float_type norm = 0.;
    Matrix<float_type> Kshape(4, 4);
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      const float_type intfac =
          weights[ig] * shape_distribution(axis_ratio[ig]);
      norm += intfac;

      const TMatrix *Tmatrix = TMatrixCalculator::calculate_TMatrix(
          1., axis_ratio[ig], 2.e-7, 1.e-4, 200, 1.e-4, 2,
          std::complex<float_type>(4., 0.1), 500);
      const Matrix<float_type> Ksingle =
          Tmatrix->get_extinction_matrix(0., 0., 0.3 * M_PI, 0.);
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          Kshape(i, j) += intfac * Ksingle(i, j);
        }
      }
      delete Tmatrix;
    }
    const float_type norm_inv = 1. / norm;
    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        Kshape(i, j) *= norm_inv;
      }
    }

    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        ctm_warning("Kshape(%" PRIuFAST8 ",%" PRIuFAST8 "): %g", i, j,
                    double(Kshape(i, j)));
      }
    }
  }

  /// test to verify the fixed orientation get_scattering_matrix() function
  if (do_fixed_orientation_test) {
    const TMatrix *Tmatrix = TMatrixCalculator::calculate_TMatrix(
        0.1, 0.5, 1.e-5, 2.e-6 * M_PI, 200, 1.e-4, 2,
        std::complex<float_type>(1.5, 0.02), 500);

    for (uint_fast32_t i = 0; i < 100; ++i) {

      const double theta_in = Utilities::random_double() * M_PI;
      const double phi_in = Utilities::random_double() * 2. * M_PI;
      const double theta_out = Utilities::random_double() * M_PI;
      const double phi_out = Utilities::random_double() * 2. * M_PI;

      Matrix<std::complex<float_type>> Sslow =
          Tmatrix->get_forward_scattering_matrix(0., 0., theta_in, phi_in,
                                                 theta_out, phi_out);
      Matrix<std::complex<float_type>> Sfast =
          Tmatrix->get_forward_scattering_matrix(theta_in, phi_in, theta_out,
                                                 phi_out);

      assert_values_equal_rel(double(Sslow(0, 0).real()),
                              double(Sfast(0, 0).real()), 1.e-6);
      assert_values_equal_rel(double(Sslow(0, 0).imag()),
                              double(Sfast(0, 0).imag()), 1.e-6);
      assert_values_equal_rel(double(Sslow(0, 1).real()),
                              double(Sfast(0, 1).real()), 1.e-6);
      assert_values_equal_rel(double(Sslow(0, 1).imag()),
                              double(Sfast(0, 1).imag()), 1.e-6);
      assert_values_equal_rel(double(Sslow(1, 0).real()),
                              double(Sfast(1, 0).real()), 1.e-6);
      assert_values_equal_rel(double(Sslow(1, 0).imag()),
                              double(Sfast(1, 0).imag()), 1.e-6);
      assert_values_equal_rel(double(Sslow(1, 1).real()),
                              double(Sfast(1, 1).real()), 1.e-6);
      assert_values_equal_rel(double(Sslow(1, 1).imag()),
                              double(Sfast(1, 1).imag()), 1.e-6);
    }
  }

  return 0;
}
