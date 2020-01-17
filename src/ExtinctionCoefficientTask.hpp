/**
 * @file ExtinctionCoefficientTask.hpp
 *
 * @brief Task that computes ExtinctionCoefficients.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef EXTINCTIONCOEFFICIENTTASK_HPP
#define EXTINCTIONCOEFFICIENTTASK_HPP

#include "Configuration.hpp"
#include "InteractionResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "Result.hpp"
#include "TMatrixResource.hpp"

#include <vector>

/**
 * @brief Result of an extinction coefficient calculation.
 */
class ExtinctionCoefficientResult : public Result {

  /*! @brief Give access to the computation task. */
  friend class ExtinctionCoefficientTask;

  /*! @brief Give access to the averaging task. */
  friend class ExtinctionShapeAveragingTask;

private:
  /*! @brief Extinction coefficients. */
  std::vector<float_type> _Qext;

  /*! @brief Polarised extinction coefficients. */
  std::vector<float_type> _Qextpol;

  /*! @brief Circularly polarised extinction coefficients. */
  std::vector<float_type> _Qextcpol;

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
  inline ExtinctionCoefficientResult(const int_fast32_t composition,
                                     const float_type size,
                                     const float_type wavelength,
                                     const uint_fast32_t number_of_angles)
      : Result(composition, size, wavelength,
               RESULTTYPE_EXTINCTIONCOEFFICIENTS),
        _Qext(number_of_angles, 0.), _Qextpol(number_of_angles, 0.),
        _Qextcpol(number_of_angles, 0.) {}

  virtual ~ExtinctionCoefficientResult() {}

  /**
   * @brief Get the size in memory of a hypothetical ExtinctionCoefficientResult
   * object with the given parameters.
   *
   * @param number_of_angles Number of angular points at which the coefficients
   * are computed.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t number_of_angles) {
    size_t size = sizeof(ExtinctionCoefficientResult);
    size += 3 * number_of_angles * sizeof(float_type);
    return size;
  }

  /**
   * @brief Get the extinction coefficient for the given index.
   *
   * @param index Index.
   * @return Corresponding extinction coefficient.
   */
  inline float_type get_Qext(const uint_fast32_t index) const {
    ctm_assert(index < _Qext.size());
    return _Qext[index];
  }

  /**
   * @brief Get the polarised extinction coefficient for the given index.
   *
   * @param index Index.
   * @return Corresponding polarised extinction coefficient.
   */
  inline float_type get_Qextpol(const uint_fast32_t index) const {
    ctm_assert(index < _Qextpol.size());
    return _Qextpol[index];
  }

  /**
   * @brief Get the circularly polarised extinction coefficient for the given
   * index.
   *
   * @param index Index.
   * @return Corresponding circularly polarised extinction coefficient.
   */
  inline float_type get_Qextcpol(const uint_fast32_t index) const {
    ctm_assert(index < _Qextcpol.size());
    return _Qextcpol[index];
  }
};

/**
 * @brief Angular grid used to compute extinction cross sections.
 */
class ExtinctionCoefficientGrid : public Resource,
                                  public Task,
                                  public Computable {

  /*! @brief Give access to the computation task. */
  friend class ExtinctionCoefficientTask;

  /*! @brief Give access to special Wigner D resources. */
  friend class ExtinctionSpecialWignerDResources;

private:
  /*! @brief Input zenith angles (in radians). */
  std::vector<float_type> _theta_in;

  /*! @brief Cosines of the input zenith angles. */
  std::vector<float_type> _cos_theta_in;

  /*! @brief Sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in;

  /*! @brief Inverse sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in_inverse;

public:
  /**
   * @brief Constructor.
   *
   * @param ntheta Number of zenith angles.
   * @param theta Grid of zenith angles (in radians, of size ntheta or more).
   */
  inline ExtinctionCoefficientGrid(const uint_fast32_t ntheta,
                                   const float_type *theta)
      : _theta_in(ntheta), _cos_theta_in(ntheta), _sin_theta_in(ntheta),
        _sin_theta_in_inverse(ntheta) {

    for (uint_fast32_t i = 0; i < ntheta; ++i) {
      _theta_in[i] = theta[i];
    }
  }

  virtual ~ExtinctionCoefficientGrid() {}

  /**
   * @brief Get the size in memory of a hypothetical ExtinctionCoefficientGrid
   * object with the given parameters.
   *
   * @param ntheta Number of zenith angles.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t ntheta) {
    size_t size = sizeof(ExtinctionCoefficientGrid);
    size += 4 * ntheta * sizeof(float_type);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _cos_theta_in.size();
    for (uint_fast32_t igauss = 0; igauss < ntheta; ++igauss) {
      _cos_theta_in[igauss] = cos(_theta_in[igauss]);
      _sin_theta_in[igauss] =
          sqrt((1. - _cos_theta_in[igauss]) * (1. + _cos_theta_in[igauss]));
      if (_sin_theta_in[igauss] != 0.) {
        _sin_theta_in_inverse[igauss] = 1. / _sin_theta_in[igauss];
      } else {
        _sin_theta_in_inverse[igauss] = 9000.;
      }
    }
  }

  /**
   * @brief Get the number of angles in the grid.
   *
   * @return Number of angles.
   */
  inline uint_fast32_t get_number_of_angles() const { return _theta_in.size(); }
};

/**
 * @brief Precomputed special Wigner D functions that depend on a specific value
 * of @f$n_{max}@f$ and @f$n_{GL}@f$.
 */
class ExtinctionSpecialWignerDResources : public Resource,
                                          public Task,
                                          public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Extinction coefficient grid to use. */
  const ExtinctionCoefficientGrid &_grid;

  /*! @brief Wigner D functions divided by sine for input angles. */
  std::vector<Matrix<float_type>> _wigner_d_sinx;

  /*! @brief Derivatives of the Wigner D functions for input angles. */
  std::vector<Matrix<float_type>> _dwigner_d;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param grid Extinction coefficient grid.
   */
  inline ExtinctionSpecialWignerDResources(
      const uint_fast32_t nmax, const ExtinctionCoefficientGrid &grid)
      : _nmax(nmax), _grid(grid) {

    const uint_fast32_t np1 = nmax + 1;

    _wigner_d_sinx.reserve(np1 * nmax);
    for (uint_fast32_t i = 0; i < np1 * nmax; ++i) {
      _wigner_d_sinx.push_back(
          Matrix<float_type>(grid._cos_theta_in.size(), nmax));
    }
    _dwigner_d.reserve(np1 * nmax);
    for (uint_fast32_t i = 0; i < np1 * nmax; ++i) {
      _dwigner_d.push_back(Matrix<float_type>(grid._cos_theta_in.size(), nmax));
    }
  }

  virtual ~ExtinctionSpecialWignerDResources() {}

  /**
   * @brief Get the size in memory of a hypothetical
   * ExtinctionSpecialWignerDResources object with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param grid Extinction coefficient grid.
   * @return Size in bytes of the hypothetical object.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const ExtinctionCoefficientGrid &grid) {
    size_t size = sizeof(ExtinctionSpecialWignerDResources);
    // input angles
    size += 2 * (nmax + 1) * nmax * grid._cos_theta_in.size() * nmax *
            sizeof(float_type);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);

    // read access
    quicksched.link_task_and_resource(*this, _grid, false);
  }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    const uint_fast32_t ntheta_in = _grid._cos_theta_in.size();
    for (uint_fast32_t m = 0; m < _nmax + 1; ++m) {
      for (uint_fast32_t n = 0; n < _nmax; ++n) {
        for (uint_fast32_t itheta_in = 0; itheta_in < ntheta_in; ++itheta_in) {
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_in[itheta_in], _grid._sin_theta_in[itheta_in],
              _grid._sin_theta_in_inverse[itheta_in], n + 1, m,
              &_wigner_d_sinx[m * _nmax + n].get_row(itheta_in)[0],
              &_dwigner_d[m * _nmax + n].get_row(itheta_in)[0]);
        }
      }
    }
    make_available();
  }

  /**
   * @brief Get the special Wigner D function for the given input angle.
   *
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding special Wigner D function value.
   */
  inline float_type get_wigner_d_sinx(const uint_fast32_t m,
                                      const uint_fast32_t itheta_in,
                                      const uint_fast32_t n,
                                      const uint_fast32_t nmax) const {

    ctm_assert(m < _wigner_d_sinx.size());
    ctm_assert(itheta_in <
               _wigner_d_sinx[m * _nmax + nmax].get_number_of_rows());
    ctm_assert(n > 0);
    ctm_assert(n - 1 <
               _wigner_d_sinx[m * _nmax + nmax].get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _wigner_d_sinx[m * _nmax + nmax](itheta_in, n - 1);
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given input
   * angle.
   *
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast32_t m,
                                  const uint_fast32_t itheta_in,
                                  const uint_fast32_t n,
                                  const uint_fast32_t nmax) const {

    ctm_assert(m < _dwigner_d.size());
    ctm_assert(itheta_in < _dwigner_d[m * _nmax + nmax].get_number_of_rows());
    ctm_assert(n > 0);
    ctm_assert(n - 1 < _dwigner_d[m * _nmax + nmax].get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _dwigner_d[m * _nmax + nmax](itheta_in, n - 1);
  }
};

/**
 * @brief Task that computes ExtinctionCoefficients.
 */
class ExtinctionCoefficientTask : public Task {
private:
  /*! @brief Zenith angle grid to use. */
  const ExtinctionCoefficientGrid &_grid;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief T-matrix to use (read only). */
  const TMatrixResource &_Tmatrix;

  /*! @brief N based resources to use (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Special Wigner D resources to use (read only). */
  const ExtinctionSpecialWignerDResources &_wigner_d;

  /*! @brief Resource in which the result is stored. */
  ExtinctionCoefficientResult &_result;

public:
  /**
   * @brief Constructor.
   *
   * @param grid Zenith angle grid to use.
   * @param interaction_variables Interaction variables.
   * @param Tmatrix T-matrix to use.
   * @param nfactors N based resources to use.
   * @param wigner_d Special Wigner D resources to use.
   * @param result Resource in which the result is stored.
   */
  inline ExtinctionCoefficientTask(
      const ExtinctionCoefficientGrid &grid,
      const InteractionVariables &interaction_variables,
      const TMatrixResource &Tmatrix, const NBasedResources &nfactors,
      const ExtinctionSpecialWignerDResources &wigner_d,
      ExtinctionCoefficientResult &result)
      : _grid(grid), _interaction_variables(interaction_variables),
        _Tmatrix(Tmatrix), _nfactors(nfactors), _wigner_d(wigner_d),
        _result(result) {}

  virtual ~ExtinctionCoefficientTask() {}

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _result, true);

    // read access
    quicksched.link_task_and_resource(*this, _interaction_variables, false);
    quicksched.link_task_and_resource(*this, _Tmatrix, false);
    quicksched.link_task_and_resource(*this, _nfactors, false);
    quicksched.link_task_and_resource(*this, _grid, false);
    quicksched.link_task_and_resource(*this, _wigner_d, false);
  }

  /**
   * @brief Get the forward scattering matrix @f$S@f$ for a scattering event
   * from the given input angles with its symmetry axis fixed to the
   * @f$z@f$-axis of the reference frame.
   *
   * @param itheta_in Index of the input zenith angle.
   * @param cosphi_out Cosine of the output azimuth angle,
   * @f$\cos(\phi{}_s)@f$.
   * @param sinphi_out Sine of the output azimuth angle,
   * @f$\sin(\phi{}_s)@f$.
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<std::complex<float_type>>
  get_forward_scattering_matrix(const uint_fast32_t itheta_in,
                                const float_type cosphi_out,
                                const float_type sinphi_out) const {

    const uint_fast32_t nmax = _Tmatrix.get_nmax();

    // now compute the matrix S^P
    // we precompute e^{i(phi_out-phi_in)}
    const std::complex<float_type> expiphi_p_out_m_in(cosphi_out, sinphi_out);
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

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      // now perform the actual sums over n and n'
      for (uint_fast32_t nn = nmin; nn < nmax + 1; ++nn) {

        // get the specific pi and tau for this n'
        const float_type pi_nn =
            m * _wigner_d.get_wigner_d_sinx(m, itheta_in, nn, nmax - 1);
        const float_type tau_nn =
            _wigner_d.get_dwigner_d(m, itheta_in, nn, nmax - 1);

        for (uint_fast32_t n = nmin; n < nmax + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n =
              m * _wigner_d.get_wigner_d_sinx(m, itheta_in, n, nmax - 1);
          const float_type tau_n =
              _wigner_d.get_dwigner_d(m, itheta_in, n, nmax - 1);

          // get the c factor for these values of n and n'
          // note the order of the arguments (took an awful long time to spot
          // a bug caused by changing the order)
          const std::complex<float_type> c_nnn = _nfactors.get_cnn(nn, n);

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
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _grid._theta_in.size();

    for (uint_fast32_t itheta_in = 0; itheta_in < ntheta; ++itheta_in) {

      Matrix<std::complex<float_type>> S =
          get_forward_scattering_matrix(itheta_in, 1., 0.);

      const float_type prefactor =
          2. * M_PI / _interaction_variables.get_wavenumber();
      _result._Qext[itheta_in] = prefactor * (S(0, 0) + S(1, 1)).imag();
      _result._Qextpol[itheta_in] = prefactor * (S(0, 0) - S(1, 1)).imag();
      _result._Qextcpol[itheta_in] = prefactor * (S(1, 1) - S(0, 0)).real();

      // normalise the coefficients
      const float_type a = _interaction_variables.get_equal_volume_radius();
      const float_type norm = M_PI * a * a;
      _result._Qext[itheta_in] /= norm;
      _result._Qextpol[itheta_in] /= norm;
      _result._Qextcpol[itheta_in] /= norm;
    }
  }
};

#endif // EXTINCTIONCOEFFICIENTTASK_HPP
