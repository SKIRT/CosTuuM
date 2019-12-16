/**
 * @file AbsorptionCoefficientTask.hpp
 *
 * @brief Task that computes AbsorptionCoefficients.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ABSORPTIONCOEFFICIENTTASK_HPP
#define ABSORPTIONCOEFFICIENTTASK_HPP

#include "InteractionResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "Result.hpp"
#include "TMatrixResource.hpp"

#include <vector>

/**
 * @brief Result of an absorption coefficient calculation.
 */
class AbsorptionCoefficientResult : public Result {

  /*! @brief Give access to the computation task. */
  friend class AbsorptionCoefficientTask;

private:
  /*! @brief Absorption coefficients. */
  std::vector<float_type> _Qabs;

  /*! @brief Polarised absorption coefficients. */
  std::vector<float_type> _Qabspol;

public:
  /**
   * @brief Constructor.
   *
   * @param composition Composition parameter value for the result.
   * @param size Particle size parameter value for the result (in m).
   * @param wavelength Wavelength value for the result (in m).
   * @param number_of_angles Number of angular points at which the coefficients
   * are computed.
   */
  inline AbsorptionCoefficientResult(const int_fast32_t composition,
                                     const float_type size,
                                     const float_type wavelength,
                                     const uint_fast32_t number_of_angles)
      : Result(composition, size, wavelength,
               RESULTTYPE_ABSORPTIONCOEFFICIENTS),
        _Qabs(number_of_angles, 0.), _Qabspol(number_of_angles, 0.) {}

  virtual ~AbsorptionCoefficientResult() {}

  /**
   * @brief Get the size in memory of a hypothetical AbsorptionCoefficientResult
   * object with the given parameters.
   *
   * @param number_of_angles Number of angular points at which the coefficients
   * are computed.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t number_of_angles) {
    size_t size = sizeof(AbsorptionCoefficientResult);
    size += 2 * number_of_angles * sizeof(float_type);
    return size;
  }

  /**
   * @brief Get the absorption coefficient for the given index.
   *
   * @param index Index.
   * @return Corresponding absorption coefficient.
   */
  inline float_type get_Qabs(const uint_fast32_t index) const {
    ctm_assert(index < _Qabs.size());
    return _Qabs[index];
  }

  /**
   * @brief Get the polarised absorption coefficient for the given index.
   *
   * @param index Index.
   * @return Corresponding polarised absorption coefficient.
   */
  inline float_type get_Qabspol(const uint_fast32_t index) const {
    ctm_assert(index < _Qabspol.size());
    return _Qabspol[index];
  }
};

/**
 * @brief Angular grid used to compute absorption cross sections.
 */
class AbsorptionCoefficientGrid {

  /*! @brief Give access to the computation task. */
  friend class AbsorptionCoefficientTask;

private:
  /*! @brief Zenith angles (in radians). */
  std::vector<float_type> _theta;

public:
  /**
   * @brief Constructor.
   *
   * @param ntheta Number of zenith angles.
   * @param theta Grid of zenith angles (in radians, of size ntheta or more).
   */
  inline AbsorptionCoefficientGrid(const uint_fast32_t ntheta,
                                   const float_type *theta)
      : _theta(ntheta) {

    for (uint_fast32_t i = 0; i < ntheta; ++i) {
      _theta[i] = theta[i];
    }
  }

  /**
   * @brief Get the number of angles in the grid.
   *
   * @return Number of angles.
   */
  inline uint_fast32_t get_number_of_angles() const { return _theta.size(); }
};

/**
 * @brief Task that computes AbsorptionCoefficients.
 */
class AbsorptionCoefficientTask : public Task {
private:
  /*! @brief Number of Gauss-Legendre quadrature points used to compute angular
   *  averages. */
  const uint_fast32_t _ngauss;

  /*! @brief Zenith angle grid to use. */
  const AbsorptionCoefficientGrid &_grid;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief T-matrix to use (read only). */
  const TMatrixResource &_Tmatrix;

  /*! @brief Resource in which the result is stored. */
  AbsorptionCoefficientResult &_result;

public:
  /**
   * @brief Constructor.
   *
   * @param ngauss Number of Gauss-Legendre quadrature points used to compute
   * angular averages.
   * @param grid Zenith angle grid to use.
   * @param interaction_variables Interaction variables.
   * @param Tmatrix T-matrix to use.
   * @param result Resource in which the result is stored.
   */
  inline AbsorptionCoefficientTask(
      const uint_fast32_t ngauss, const AbsorptionCoefficientGrid &grid,
      const InteractionVariables &interaction_variables,
      const TMatrixResource &Tmatrix, AbsorptionCoefficientResult &result)
      : _ngauss(ngauss), _grid(grid),
        _interaction_variables(interaction_variables), _Tmatrix(Tmatrix),
        _result(result) {}

  virtual ~AbsorptionCoefficientTask() {}

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _result, true);
  }

  /**
   * @brief Get the forward scattering matrix @f$S@f$ for a scattering event
   * from the given input angles to the given output angles at a particle with
   * its symmetry axis fixed to the @f$z@f$-axis of the reference frame.
   *
   * @param costheta_in Cosine of the input zenith angle,
   * @f$\cos(\theta{}_i)@f$.
   * @param sintheta_in Sine of the input zenith angle,
   * @f$\sin(\theta{}_i)@f$.
   * @param sintheta_in_inv Inverse sine of the input zenith angle,
   * @f$\frac{1}{\sin(\theta{}_i)}@f$.
   * @param cosphi_in Cosine of the input azimuth angle,
   * @f$\cos(\phi{}_i)@f$.
   * @param sinphi_in Sine of the input azimuth angle,
   * @f$\sin(\phi{}_i)@f$.
   * @param costheta_out Cosine of the output zenith angle,
   * @f$\cos(\theta{}_s)@f$.
   * @param sintheta_out Sine of the output zenith angle,
   * @f$\sin(\theta{}_s)@f$.
   * @param sintheta_out_inv Inverse sine of the output zenith angle,
   * @f$\frac{1}{\sin(\theta{}_s)}@f$.
   * @param cosphi_out Cosine of the output azimuth angle,
   * @f$\cos(\phi{}_s)@f$.
   * @param sinphi_out Sine of the output azimuth angle,
   * @f$\sin(\phi{}_s)@f$.
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<std::complex<float_type>> get_forward_scattering_matrix(
      const float_type costheta_in, const float_type sintheta_in,
      const float_type sintheta_in_inv, const float_type cosphi_in,
      const float_type sinphi_in, const float_type costheta_out,
      const float_type sintheta_out, const float_type sintheta_out_inv,
      const float_type cosphi_out, const float_type sinphi_out) const {

    const uint_fast32_t nmax = _Tmatrix.get_nmax();
    // precompute the c factors
    const std::complex<float_type> icompl(0., 1.);
    Matrix<std::complex<float_type>> c(nmax, nmax);
    std::complex<float_type> icomp_pow_nn = icompl;
    for (uint_fast32_t nn = 1; nn < nmax + 1; ++nn) {
      std::complex<float_type> icomp_pow_m_n_m_1(-1.);
      for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
        // icomp_pow_nn*icomp_pow_m_n_m_1 now equals i^(nn - n - 1)
        c(n - 1, nn - 1) = icomp_pow_m_n_m_1 * icomp_pow_nn *
                           float_type(sqrt((2. * n + 1.) * (2. * nn + 1.) /
                                           (n * nn * (n + 1.) * (nn + 1.))));
        icomp_pow_m_n_m_1 /= icompl;
      }
      icomp_pow_nn *= icompl;
    }

    // now compute the matrix S^P
    // we precompute e^{i(phi_out-phi_in)}
    const std::complex<float_type> expiphi_p_out_m_in(
        cosphi_out * cosphi_in + sinphi_out * sinphi_in,
        sinphi_out * cosphi_in - cosphi_out * sinphi_in);
    // e^{im(phi_out-phi_in)} is computed recursively, starting with the value
    // for m=0: 1
    std::complex<float_type> expimphi_p_out_m_in(1., 0.);
    Matrix<std::complex<float_type>> S(2, 2);
    // instead of summing over n and n', we sum over m, since then we can reuse
    // the e^{im(phi_out-phi_in)}, pi and tau factors
    for (uint_fast32_t m = 0; m < nmax + 1; ++m) {
      // only n and n' values larger than or equal to m have non-trivial
      // contributions to the S matrix
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));

      // precompute the pi and tau functions for this value of m
      std::vector<float_type> pi_in(nmax), tau_in(nmax);
      SpecialFunctions::wigner_dn_0m_sinx(costheta_in, sintheta_in,
                                          sintheta_in_inv, nmax, m, &pi_in[0],
                                          &tau_in[0]);
      std::vector<float_type> pi_out(nmax), tau_out(nmax);
      SpecialFunctions::wigner_dn_0m_sinx(costheta_out, sintheta_out,
                                          sintheta_out_inv, nmax, m, &pi_out[0],
                                          &tau_out[0]);

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      // now perform the actual sums over n and n'
      for (uint_fast32_t nn = nmin; nn < nmax + 1; ++nn) {

        // get the specific pi and tau for this n'
        const float_type pi_nn = m * pi_in[nn - 1];
        const float_type tau_nn = tau_in[nn - 1];

        for (uint_fast32_t n = nmin; n < nmax + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n = m * pi_out[n - 1];
          const float_type tau_n = tau_out[n - 1];

          // get the c factor for these values of n and n'
          const std::complex<float_type> c_nnn = c(n - 1, nn - 1);

          // get the T11 and T22 elements for this m, n and n' (we need these
          // in all cases)
          const std::complex<float_type> T11nmnnm = _Tmatrix(0, n, m, 0, nn, m);
          const std::complex<float_type> T22nmnnm = _Tmatrix(1, n, m, 1, nn, m);
          // if m=0, the T12 and T21 matrices are trivially zero, and we can
          // simplify the expression for S
          if (m == 0) {
            const std::complex<float_type> factor = c_nnn * tau_n * tau_nn;
            S(0, 0) += factor * T22nmnnm;
            S(1, 1) += factor * T11nmnnm;
          } else {
            // in the general case m=/=0, we also need the T12 and T21 elements
            // for this m, n and n'
            const std::complex<float_type> T12nmnnm =
                _Tmatrix(0, n, m, 1, nn, m);
            const std::complex<float_type> T21nmnnm =
                _Tmatrix(1, n, m, 0, nn, m);

            // due to m symmetry, S11 and S22 only have the cosine factor,
            // while S12 and S21 only have the sine factor
            const std::complex<float_type> real_factor = c_nnn * fcos;
            const std::complex<float_type> imag_factor = c_nnn * fsin;

            // precompute the pi and tau factor combinations
            const float_type pi_pi = pi_n * pi_nn;
            const float_type pi_tau = pi_n * tau_nn;
            const float_type tau_pi = tau_n * pi_nn;
            const float_type tau_tau = tau_n * tau_nn;

            S(0, 0) += real_factor * (T11nmnnm * pi_pi + T21nmnnm * tau_pi +
                                      T12nmnnm * pi_tau + T22nmnnm * tau_tau);
            S(0, 1) += imag_factor * (T11nmnnm * pi_tau + T21nmnnm * tau_tau +
                                      T12nmnnm * pi_pi + T22nmnnm * tau_pi);
            S(1, 0) -= imag_factor * (T11nmnnm * tau_pi + T21nmnnm * pi_pi +
                                      T12nmnnm * tau_tau + T22nmnnm * pi_tau);
            S(1, 1) += real_factor * (T11nmnnm * tau_tau + T21nmnnm * pi_tau +
                                      T12nmnnm * tau_pi + T22nmnnm * pi_pi);
          }
        }
      }
    }
    // now divide all expressions by the wavenumber
    const float_type kinv = 1. / _interaction_variables.get_wavenumber();
    S(0, 0) *= kinv;
    S(0, 1) *= kinv;
    S(1, 0) *= kinv;
    S(1, 1) *= kinv;

    return S;
  }

  /**
   * @brief Get the forward scattering matrix @f$S@f$ for a scattering event
   * from the given input angles to the given output angles at a particle with
   * its symmetry axis fixed to the @f$z@f$-axis of the reference frame.
   *
   * @param theta_in_radians Zenith angle of the incoming photon,
   * @f$\theta{}_i@f$ (in radians).
   * @param phi_in_radians Azimuth angle of the incoming photon, @f$\phi{}_i@f$
   * (in radians).
   * @param theta_out_radians Zenith angle of the scattered photon,
   * @f$\theta{}_s@f$ (in radians).
   * @param phi_out_radians Azimuth angle fo the scattered photon,
   * @f$\phi{}_s@f$ (in radians).
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<std::complex<float_type>>
  get_forward_scattering_matrix(const float_type theta_in_radians,
                                const float_type phi_in_radians,
                                const float_type theta_out_radians,
                                const float_type phi_out_radians) const {

    // Mishchenko includes some (buggy) corrections for small angles
    // might be worth looking into this in a later stage...

    // compute all sines and cosines in one go; we need all of them anyway
    const float_type costheta_l_in = cos(theta_in_radians);
    const float_type sintheta_l_in = sin(theta_in_radians);
    const float_type costheta_l_out = cos(theta_out_radians);
    const float_type sintheta_l_out = sin(theta_out_radians);
    const float_type cosphi_l_in = cos(phi_in_radians);
    const float_type sinphi_l_in = sin(phi_in_radians);
    const float_type cosphi_l_out = cos(phi_out_radians);
    const float_type sinphi_l_out = sin(phi_out_radians);

    float_type sintheta_l_in_inv;
    if (sintheta_l_in != 0.) {
      sintheta_l_in_inv = 1. / sintheta_l_in;
    } else {
      // this value will not be used, but we set it to something anyway
      sintheta_l_in_inv = 9000.;
    }

    float_type sintheta_l_out_inv;
    if (sintheta_l_out != 0.) {
      sintheta_l_out_inv = 1. / sintheta_l_out;
    } else {
      // this value will not be used, but we set it to something anyway
      sintheta_l_out_inv = 9000.;
    }

    return get_forward_scattering_matrix(
        costheta_l_in, sintheta_l_in, sintheta_l_in_inv, cosphi_l_in,
        sinphi_l_in, costheta_l_out, sintheta_l_out, sintheta_l_out_inv,
        cosphi_l_out, sinphi_l_out);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    std::vector<float_type> thetaGL(_ngauss), costhetaGL(_ngauss),
        costhetaweightsGL(_ngauss);
    SpecialFunctions::get_gauss_legendre_points_and_weights<float_type>(
        _ngauss, costhetaGL, costhetaweightsGL);
    std::vector<float_type> phiGL(_ngauss), phiweightsGL(_ngauss);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        _ngauss, 0., 2. * M_PI, phiGL, phiweightsGL);

    for (uint_fast32_t igauss = 0; igauss < _ngauss; ++igauss) {
      thetaGL[igauss] = acos(costhetaGL[igauss]);
    }

    const uint_fast32_t ntheta = _grid._theta.size();
    for (uint_fast32_t itheta_in = 0; itheta_in < ntheta; ++itheta_in) {

      const float_type theta = _grid._theta[itheta_in];
      Matrix<std::complex<float_type>> S =
          get_forward_scattering_matrix(theta, 0., theta, 0.);

      const float_type prefactor =
          2. * M_PI / _interaction_variables.get_wavenumber();
      _result._Qabs[itheta_in] = prefactor * (S(0, 0) + S(1, 1)).imag();
      _result._Qabspol[itheta_in] = prefactor * (S(0, 0) - S(1, 1)).imag();

      const float_type half(0.5);
      for (uint_fast32_t itheta = 0; itheta < _ngauss; ++itheta) {
        const float_type theta_out = thetaGL[itheta];
        for (uint_fast32_t iphi = 0; iphi < _ngauss; ++iphi) {
          const float_type phi_out = phiGL[iphi];

          Matrix<std::complex<float_type>> Stp =
              get_forward_scattering_matrix(theta, 0., theta_out, phi_out);

          const float_type weight =
              costhetaweightsGL[itheta] * phiweightsGL[iphi];
          const float_type Z00 =
              (half *
               (Stp(0, 0) * conj(Stp(0, 0)) + Stp(0, 1) * conj(Stp(0, 1)) +
                Stp(1, 0) * conj(Stp(1, 0)) + Stp(1, 1) * conj(Stp(1, 1))))
                  .real();
          _result._Qabs[itheta_in] -= Z00 * weight;

          const float_type Z10 =
              (half *
               (Stp(0, 0) * conj(Stp(0, 0)) + Stp(0, 1) * conj(Stp(0, 1)) -
                Stp(1, 0) * conj(Stp(1, 0)) - Stp(1, 1) * conj(Stp(1, 1))))
                  .real();
          _result._Qabspol[itheta_in] -= Z10 * weight;
        }
      }

      // normalise the coefficients
      const float_type a = _interaction_variables.get_equal_volume_radius();
      const float_type norm = M_PI * a * a;
      _result._Qabs[itheta_in] /= norm;
      _result._Qabspol[itheta_in] /= norm;
    }
  }
};

#endif // ABSORPTIONCOEFFICIENTTASK_HPP
