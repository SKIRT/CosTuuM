/**
 * @file TMatrixCalculator.hpp
 *
 * @brief Class that calculates the T-matrix for a single particle.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TMATRIXCALCULATOR_HPP
#define TMATRIXCALCULATOR_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "OrientationDistribution.hpp"
#include "SpecialFunctions.hpp"
#include "TMatrix.hpp"

/**
 * @brief Class that calculates the T-matrix for a single particle.
 */
class TMatrixCalculator {

public:
  /**
   * @brief Calculate the T-matrix for the particle with the given properties
   * and the given wavelength.
   *
   * @param ratio_of_radii Ratio of the input radius to the equal volume sphere
   * radius. If this value is 1, this signals that the given particle radius is
   * the equal volume sphere radius. If not, the value is assumed to be the
   * equal surface area sphere radius and this ratio is updated by the function.
   * @param axis_ratio Ratio of the horizontal and vertical axis of the
   * spheroidal particle.
   * @param axi Radius of the particle (in m).
   * @param wavelength Wavelength of the incident radiation (in m).
   * @param maximum_order Maximum allowed order of the T-matrix expansion.
   * @param tolerance Tolerance for the calculation. The expansion stops when
   * the relative difference between the scattering and extinction amplitudes
   * for two successive orders drops below this value.
   * @param ndgs Number of Gauss-Legendre quadrature points per order of
   * expansion during the first phase of the calculation.
   * @param mr Complex refractive index of the scattering particle.
   * @param maximum_ngauss Maximum number of Gauss-Legendre quadrature points.
   * @return Pointer to the T-matrix. Memory management for this pointer is
   * transferred to the caller.
   */
  inline static TMatrix *
  calculate_TMatrix(float_type ratio_of_radii, const float_type axis_ratio,
                    const float_type axi, const float_type wavelength,
                    const uint_fast32_t maximum_order,
                    const float_type tolerance, const uint_fast32_t ndgs,
                    const std::complex<float_type> mr,
                    const uint_fast32_t maximum_ngauss) {

    // make sure 'ratio_of_radii' contains the right ratio if it is not 1
    if (abs(ratio_of_radii - 1.) > 1.e-8) {
      ratio_of_radii =
          SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
              axis_ratio);
    }

    // R_V is the equivalent sphere radius
    const float_type R_V = ratio_of_radii * axi;
    // we need a reasonable initial guess for the order of the spherical
    // harmonics the below expression provides this guess
    const float_type xev = 2. * M_PI * R_V / wavelength;
    uint_fast32_t nmax = static_cast<uint_fast32_t>(
        std::max(float_type(4.), xev + 4.05 * cbrt(xev)));

    // we need to find a maximum expansion order and number of Gauss-Legendre
    // quadrature points that leads to a converged T-matrix
    // To this end, we will start with a reasonable guess for the order, and
    // then keep track of how the accuracy of the scattering and extinction
    // factors changes as we first increase the order and then the number of
    // quadrature points To simplify the calculation, we only compute the m=0
    // degree terms in the T-matrix during the convergence loops
    TMatrix *active_Tmatrix = nullptr;

    // loop control variables
    float_type old_qext = 0.;
    float_type old_qsca = 0.;
    float_type dext = 1.;
    float_type dsca = 1.;
    while (nmax < maximum_order && (dext > tolerance || dsca > tolerance)) {
      // initially we assume that the number of quadrature points is a fixed
      // multiple of the order
      const uint_fast32_t ngauss = ndgs * nmax;

      // delete the old matrix (if it exists) and create a new one
      delete active_Tmatrix;
      active_Tmatrix =
          new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
      TMatrix &T = *active_Tmatrix;

      // calculate the scattering and extinction factors for this iteration
      float_type qsca = 0.;
      float_type qext = 0.;
      for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
        const float_type dn1 = 2. * n + 1.;
        qsca += dn1 * (std::norm(T(0, n, 0, 0, n, 0)) +
                       std::norm(T(1, n, 0, 1, n, 0)));
        qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
      }
      // compute the relative difference w.r.t. the previous values
      dsca = abs((old_qsca - qsca) / qsca);
      dext = abs((old_qext - qext) / qext);
      old_qext = qext;
      old_qsca = qsca;

      // some (temporary) diagnostic output
      ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

      // increase the order
      ++nmax;
    }

    // check if we managed to converge
    if (nmax == maximum_order) {
      ctm_error("Unable to converge!");
    } else {
      // correct for overshoot in last iteration
      --nmax;
      ctm_warning("Converged for nmax = %" PRIuFAST32, nmax);
    }

    // now start increasing the number of quadrature points
    uint_fast32_t ngauss = nmax * ndgs + 1;
    dext = tolerance + 1.;
    dsca = tolerance + 1.;
    while (ngauss < maximum_ngauss && (dext > tolerance || dsca > tolerance)) {

      // delete the old matrix and create a new one
      delete active_Tmatrix;
      active_Tmatrix =
          new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
      TMatrix &T = *active_Tmatrix;

      // calculate the scattering and extinction factors
      float_type qsca = 0.;
      float_type qext = 0.;
      for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
        const float_type dn1 = 2. * n + 1.;
        qsca += dn1 * (std::norm(T(0, n, 0, 0, n, 0)) +
                       std::norm(T(1, n, 0, 1, n, 0)));
        qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
      }
      // compute the relative difference w.r.t. the old values
      dsca = abs((old_qsca - qsca) / qsca);
      dext = abs((old_qext - qext) / qext);
      old_qext = qext;
      old_qsca = qsca;

      // some diagnostic output
      ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

      // increase the number of quadrature points
      ++ngauss;
    }

    // check if we reached convergence
    if (ngauss == maximum_ngauss) {
      ctm_error("Unable to converge!");
    } else {
      // correct for overshoot in final iteration
      --ngauss;
      ctm_warning("Converged for ngauss = %" PRIuFAST32, ngauss);
    }

    // ok: we found the right order and number of quadrature points
    // we now need to compute the additional elements for m=/=0
    active_Tmatrix->compute_additional_elements();

    // compute the actual scattering and extinction factors using the full
    // matrix
    TMatrix &T = *active_Tmatrix;
    old_qsca = T.get_scattering_coefficient();
    old_qext = T.get_extinction_coefficient();
    // output the factors
    ctm_warning("qsca: %g", double(old_qsca));
    ctm_warning("qext: %g", double(old_qext));
    ctm_warning("walb: %g", double(-old_qsca / old_qext));

    return active_Tmatrix;
  }

  /**
   * @brief Compute the T-matrix for an ensemble of dust particles with the
   * given individual T-matrix, distributed according to the given orientation
   * distribution.
   *
   * The method we use here is based on Mishchenko, 1991, The Astrophysical
   * Journal, 367, 561 (https://doi.org/10.1086/169652).
   *
   * @param T_single T-matrix for a single dust particle in the ensemble.
   * @param orientation_distribution Orientation distribution for the dust
   * particles in the ensemble.
   * @return Pointer to the T-matrix for the ensemble of dust particles.
   * Memory management for this pointer is transferred to the caller.
   */
  inline static TMatrix *apply_orientation_distribution(
      const TMatrix &T_single,
      const OrientationDistribution &orientation_distribution) {

    // create a copy of the original T matrix; we do not care about the old
    // values, but this way the dimensions of the ensemble T matrix will be
    // the same
    TMatrix *T_ensemble = new TMatrix(T_single);

    const uint_fast32_t nmax = T_single.get_nmax();

    // check that we precomputed enough expansion coefficients in the expansion
    // of the orientation distribution
    ctm_assert(orientation_distribution.get_maximum_order() >= 2 * nmax);

    // loop over all T matrix elements:
    //  - outer loop over m
    //  - inner loops over n1 and n2
    for (uint_fast32_t m = 0; m < nmax + 1; ++m) {
      // msign = (-1)^m
      const int_fast8_t msign = (m % 2 == 0) ? 1 : -1;
      // the minimum n value depends on the value of m
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));
      // inner loops
      for (uint_fast32_t n1 = nmin; n1 < nmax + 1; ++n1) {
        for (uint_fast32_t n2 = nmin; n2 < nmax + 1; ++n2) {
          // compute a new value for the T matrix coefficients T_{mn_1mn_2}
          std::complex<float_type> Tn1n2[2][2];
          // std::abs does not seem to work well in this context when quad
          // precision is activated, so we just mimic it ourselves
          const uint_fast32_t n12min = (n1 > n2) ? n1 - n2 : n2 - n1;
          const uint_fast32_t n12max = n1 + n2;
          const uint_fast32_t M = std::min(n1, n2);
          // first summation in Mishchenko (1991), eq. 3.27
          for (uint_fast32_t N = n12min; N < n12max + 1; ++N) {
            // now that we have n1, n2 and N, we can compute the Clebsch-Gordan
            // coefficients for this element of the sum
            const std::vector<float_type> CGcoeff =
                SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(
                    n1, n2, N);
            // we can also get the expansion coefficient
            const float_type pN = orientation_distribution.get_coefficient(N);
            // second summation in Mishchenko (1991), eq. 3.27
            for (uint_fast32_t m1 = 0; m1 < M + 1; ++m1) {
              // all elements m1=/=0 contribute twice because of symmetry
              float_type factor(2.);
              if (m1 == 0) {
                factor = 1.;
              }
              // m1sign = (-1)^m1
              const int_fast8_t m1sign = (m1 % 2 == 0) ? 1 : -1;
              const int_fast8_t mm1sign = m1sign * msign;
              ctm_assert(M + m < CGcoeff.size());
              ctm_assert(M + m1 < CGcoeff.size());
              // prefactor for all terms in this part of the sum
              const float_type CGfac =
                  mm1sign * factor * pN * CGcoeff[M + m] * CGcoeff[M + m1];
              // depending on the sign of N + n1 + n2, odd/even terms in the
              // sum cancel out because of symmetry
              if ((N + n12max) % 2 == 0) {
                Tn1n2[0][0] += CGfac * T_single(0, n1, m1, 0, n2, m1);
                Tn1n2[1][1] += CGfac * T_single(1, n1, m1, 1, n2, m1);
              } else {
                Tn1n2[0][1] += CGfac * T_single(0, n1, m1, 1, n2, m1);
                Tn1n2[1][0] += CGfac * T_single(1, n1, m1, 0, n2, m1);
              }
            }
          }
          // now set the new values for T_{mn_1mn_2}
          T_ensemble->set_element(0, n1, m, 0, n2, m, Tn1n2[0][0]);
          T_ensemble->set_element(0, n1, m, 1, n2, m, Tn1n2[0][1]);
          T_ensemble->set_element(1, n1, m, 0, n2, m, Tn1n2[1][0]);
          T_ensemble->set_element(1, n1, m, 1, n2, m, Tn1n2[1][1]);
        }
      }
    }

    return T_ensemble;
  }
};

#endif // TMATRIXCALCULATOR_HPP
