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
};

#endif // TMATRIXCALCULATOR_HPP