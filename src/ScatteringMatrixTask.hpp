/**
 * @file ScatteringMatrixTask.hpp
 *
 * @brief Task that computes a scattering matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SCATTERINGMATRIXTASK_HPP
#define SCATTERINGMATRIXTASK_HPP

#include "Configuration.hpp"
#include "InteractionResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "Result.hpp"
#include "TMatrixResource.hpp"

#include <vector>

/**
 * @brief Result of scattering matrix calculation.
 */
class ScatteringMatrixResult : public Result {

  /*! @brief Give access to the computation task. */
  friend class ScatteringMatrixTask;

  /*! @brief Give access to the averaging task. */
  friend class ScatteringMatrixShapeAveragingTask;

private:
  /*! @brief Scattering matrices. */
  std::vector<Matrix<float_type>> _scattering_matrix;

public:
  /**
   * @brief Constructor.
   *
   * @param composition Composition parameter value for the result.
   * @param size Particle size parameter value for the result (in m).
   * @param wavelength Wavelength value for the result (in m).
   * @param number_of_angles Number of angular points at which the matrix is
   * computed.
   */
  inline ScatteringMatrixResult(const int_fast32_t composition,
                                const float_type size,
                                const float_type wavelength,
                                const uint_fast32_t number_of_angles)
      : Result(composition, size, wavelength,
               RESULTTYPE_ABSORPTIONCOEFFICIENTS) {

    _scattering_matrix.reserve(number_of_angles);
    for (uint_fast32_t i = 0; i < number_of_angles; ++i) {
      _scattering_matrix.push_back(Matrix<float_type>(4, 4));
    }
  }

  virtual ~ScatteringMatrixResult() {}

  /**
   * @brief Get the size in memory of a hypothetical ScatteringMatrixResult
   * object with the given parameters.
   *
   * @param number_of_angles Number of angular points at which the matrix is
   * computed.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t number_of_angles) {
    size_t size = sizeof(ScatteringMatrixResult);
    size += 16 * number_of_angles * sizeof(float_type);
    return size;
  }

  /**
   * @brief Get the scattering matrix for the given index.
   *
   * @param index Index.
   * @return Corresponding scattering matrix.
   */
  inline const Matrix<float_type> &
  get_scattering_matrix(const uint_fast32_t index) const {
    ctm_assert(index < _scattering_matrix.size());
    return _scattering_matrix[index];
  }
};

/**
 * @brief Angular grid used to compute scattering matrices.
 */
class ScatteringMatrixGrid : public Resource, public Task, public Computable {

  /*! @brief Give access to the computation task. */
  friend class ScatteringMatrixTask;

  /*! @brief Give access to special Wigner D resources. */
  friend class ScatteringMatrixSpecialWignerDResources;

private:
  /*! @brief Input zenith angle (in radians). */
  const float_type _theta_in;

  /*! @brief Cosine of the input zenith angle. */
  const float_type _cos_theta_in;

  /*! @brief Sine of the input zenith angle. */
  const float_type _sin_theta_in;

  /*! @brief Inverse sine of the input zenith angle. */
  const float_type _sin_theta_in_inverse;

  /*! @brief Output zenith angles. */
  std::vector<float_type> _theta_out;

  /*! @brief Cosines of the output zenith angles. */
  std::vector<float_type> _cos_theta_out;

  /*! @brief Sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out;

  /*! @brief Inverse sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out_inverse;

  /*! @brief Output azimuth angles. */
  std::vector<float_type> _phi_out;

  /*! @brief Cosines of the output azimuth angles. */
  std::vector<float_type> _cos_phi_out;

  /*! @brief Sines of the output azimuth angles. */
  std::vector<float_type> _sin_phi_out;

  /*! @brief Sample the grid points using Gauss-Legendre quadrature points? */
  const bool _use_gauss_legendre_samples;

public:
  /**
   * @brief Constructor.
   *
   * @param theta_in Input zenith angle.
   * @param ntheta_out Number of output zenith angles.
   * @param nphi_out Number of output azimuth angles.
   * @param use_gauss_legendre_samples Sample the grid points using
   * Gauss-Legendre quadrature points?
   */
  inline ScatteringMatrixGrid(const float_type theta_in,
                              const uint_fast32_t ntheta_out,
                              const uint_fast32_t nphi_out,
                              const bool use_gauss_legendre_samples = true)
      : _theta_in(theta_in), _cos_theta_in(cos(theta_in)),
        _sin_theta_in(sqrt((1. - _cos_theta_in) * (1. + _cos_theta_in))),
        _sin_theta_in_inverse((_sin_theta_in != 0.) ? 1. / _sin_theta_in
                                                    : 9000.),
        _theta_out(ntheta_out, 0.), _cos_theta_out(ntheta_out, 0.),
        _sin_theta_out(ntheta_out, 0.), _sin_theta_out_inverse(ntheta_out, 0.),
        _phi_out(nphi_out, 0.), _cos_phi_out(nphi_out, 0.),
        _sin_phi_out(nphi_out, 0.),
        _use_gauss_legendre_samples(use_gauss_legendre_samples) {}

  virtual ~ScatteringMatrixGrid() {}

  /**
   * @brief Get the size in memory of a hypothetical ScatteringMatrixGrid
   * object with the given parameters.
   *
   * @param theta_in Input zenith angle.
   * @param ntheta_out Number of output zenith angles.
   * @param nphi_out Number of output azimuth angles.
   * @return Size of the hypothetical object (in bytes).
   */
  static inline size_t get_memory_size(const float_type theta_in,
                                       const uint_fast32_t ntheta_out,
                                       const uint_fast32_t nphi_out) {
    size_t size = sizeof(ScatteringMatrixGrid);
    size += 4 * ntheta_out * sizeof(float_type);
    size += 4 * nphi_out * sizeof(float_type);
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

    const uint_fast32_t ntheta_out = _theta_out.size();
    float_type dtheta = 0.;
    if (_use_gauss_legendre_samples) {
      std::vector<float_type> theta_weights_ignored(ntheta_out);
      SpecialFunctions::get_gauss_legendre_points_and_weights<float_type>(
          ntheta_out, _cos_theta_out, theta_weights_ignored);
    } else {
      dtheta = M_PI / ntheta_out;
    }
    for (uint_fast32_t i = 0; i < ntheta_out; ++i) {
      if (_use_gauss_legendre_samples) {
        _theta_out[i] = acos(_cos_theta_out[i]);
      } else {
        _theta_out[i] = i * dtheta;
        _cos_theta_out[i] = cos(_theta_out[i]);
      }
      _sin_theta_out[i] =
          sqrt((1. - _cos_theta_out[i]) * (1. + _cos_theta_out[i]));
      if (_sin_theta_out[i] != 0.) {
        _sin_theta_out_inverse[i] = 1. / _sin_theta_out[i];
      } else {
        _sin_theta_out_inverse[i] = 9000.;
      }
    }

    const uint_fast32_t nphi_out = _phi_out.size();
    float_type dphi = 0.;
    if (_use_gauss_legendre_samples) {
      std::vector<float_type> phi_weights_ignored(nphi_out);
      SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
          nphi_out, 0., 2. * M_PI, _phi_out, phi_weights_ignored);
    } else {
      dphi = 2. * M_PI / nphi_out;
    }
    for (uint_fast32_t i = 0; i < nphi_out; ++i) {
      if (!_use_gauss_legendre_samples) {
        _phi_out[i] = i * dphi;
      }
      _cos_phi_out[i] = cos(_phi_out[i]);
      _sin_phi_out[i] = sin(_phi_out[i]);
    }

    make_available();
  }

  /**
   * @brief Get the number of angles in the grid.
   *
   * @return Number of angles.
   */
  inline uint_fast32_t get_number_of_angles() const {
    return _theta_out.size() * _phi_out.size();
  }
};

/**
 * @brief Precomputed special Wigner D functions that depend on a specific value
 * of @f$n_{max}@f$ and @f$n_{GL}@f$.
 */
class ScatteringMatrixSpecialWignerDResources : public Resource,
                                                public Task,
                                                public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Absorption coefficient grid to use. */
  const ScatteringMatrixGrid &_grid;

  /*! @brief Wigner D functions divided by sine. */
  std::vector<Matrix<float_type>> _wigner_d_sinx[2];

  /*! @brief Derivatives of the Wigner D functions. */
  std::vector<Matrix<float_type>> _dwigner_d[2];

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param grid Absorption coefficient grid.
   */
  inline ScatteringMatrixSpecialWignerDResources(
      const uint_fast32_t nmax, const ScatteringMatrixGrid &grid)
      : _nmax(nmax), _grid(grid) {

    const uint_fast32_t number_of_elements = (2 + nmax + 1) * nmax / 2;

    _wigner_d_sinx[0].reserve(number_of_elements);
    _wigner_d_sinx[1].reserve(number_of_elements);
    _dwigner_d[0].reserve(number_of_elements);
    _dwigner_d[1].reserve(number_of_elements);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        _wigner_d_sinx[0].push_back(Matrix<float_type>(1, n));
        _dwigner_d[0].push_back(Matrix<float_type>(1, n));
        _wigner_d_sinx[1].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
        _dwigner_d[1].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
      }
    }
    ctm_assert(_wigner_d_sinx[0].size() == number_of_elements);
  }

  virtual ~ScatteringMatrixSpecialWignerDResources() {}

  /**
   * @brief Get the size in memory of a hypothetical SpecialWignerDResources
   * object with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param grid Scattering matrix grid.
   * @return Size in bytes of the hypothetical object.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const ScatteringMatrixGrid &grid) {
    size_t size = sizeof(ScatteringMatrixSpecialWignerDResources);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        size += 2 * n * sizeof(float_type);
        size += 2 * grid._cos_theta_out.size() * n * sizeof(float_type);
      }
    }
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

    const uint_fast32_t ntheta_out = _grid._cos_theta_out.size();
    for (uint_fast32_t n = 1; n < _nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        const uint_fast32_t index = (2 + n) * (n - 1) / 2 + m;

        // no need to loop, there is only one input angle
        ctm_assert(index < _wigner_d_sinx[0].size());
        ctm_assert_message(
            _wigner_d_sinx[0][index].get_number_of_columns() == n,
            "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
            ", index: %" PRIuFAST32,
            _wigner_d_sinx[0][index].get_number_of_columns(), n, m, index);
        SpecialFunctions::wigner_dn_0m_sinx(
            _grid._cos_theta_in, _grid._sin_theta_in,
            _grid._sin_theta_in_inverse, n, m,
            &_wigner_d_sinx[0][index].get_row(0)[0],
            &_dwigner_d[0][index].get_row(0)[0]);

        for (uint_fast32_t itheta_out = 0; itheta_out < ntheta_out;
             ++itheta_out) {
          ctm_assert(index < _wigner_d_sinx[1].size());
          ctm_assert_message(
              _wigner_d_sinx[1][index].get_number_of_columns() == n,
              "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
              ", index: %" PRIuFAST32,
              _wigner_d_sinx[1][index].get_number_of_columns(), n, m, index);
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_out[itheta_out],
              _grid._sin_theta_out[itheta_out],
              _grid._sin_theta_out_inverse[itheta_out], n, m,
              &_wigner_d_sinx[1][index].get_row(itheta_out)[0],
              &_dwigner_d[1][index].get_row(itheta_out)[0]);
        }
      }
    }
    make_available();
  }

  /**
   * @brief Get the special Wigner D function for the given input angle.
   *
   * @param igrid Internal grid to sample.
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding special Wigner D function value.
   */
  inline float_type get_wigner_d_sinx(const uint_fast8_t igrid,
                                      const uint_fast32_t m,
                                      const uint_fast32_t itheta_in,
                                      const uint_fast32_t n,
                                      const uint_fast32_t nmax) const {

    ctm_assert(nmax > 0);
    const uint_fast32_t index = (2 + nmax) * (nmax - 1) / 2 + m;
    ctm_assert(index < _wigner_d_sinx[igrid].size());
    ctm_assert(n > 0);
    ctm_assert(m <= n);
    ctm_assert(itheta_in < _wigner_d_sinx[igrid][index].get_number_of_rows());
    ctm_assert_message(n - 1 <
                           _wigner_d_sinx[igrid][index].get_number_of_columns(),
                       "m: %" PRIuFAST32 ", n: %" PRIuFAST32
                       ", nmax: %" PRIuFAST32 ", index: %" PRIuFAST32,
                       m, n, nmax, index);
    // check that the resource was actually computed
    check_use();
    return _wigner_d_sinx[igrid][index](itheta_in, n - 1);
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given input
   * angle.
   *
   * @param igrid Internal grid to sample.
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast8_t igrid,
                                  const uint_fast32_t m,
                                  const uint_fast32_t itheta_in,
                                  const uint_fast32_t n,
                                  const uint_fast32_t nmax) const {

    ctm_assert(nmax > 0);
    const uint_fast32_t index = (2 + nmax) * (nmax - 1) / 2 + m;
    ctm_assert(index < _dwigner_d[igrid].size());
    ctm_assert(n > 0);
    ctm_assert(m <= n);
    ctm_assert(itheta_in < _dwigner_d[igrid][index].get_number_of_rows());
    ctm_assert(n - 1 < _dwigner_d[igrid][index].get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _dwigner_d[igrid][index](itheta_in, n - 1);
  }
};

/**
 * @brief Task that computes the scattering matrix.
 */
class ScatteringMatrixTask : public Task {
private:
  /*! @brief Angular grid to use. */
  const ScatteringMatrixGrid &_grid;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief T-matrix to use (read only). */
  const TMatrixResource &_Tmatrix;

  /*! @brief N based resources to use (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Special Wigner D resources to use (read only). */
  const ScatteringMatrixSpecialWignerDResources &_wigner_d;

  /*! @brief Resource in which the result is stored. */
  ScatteringMatrixResult &_result;

public:
  /**
   * @brief Constructor.
   *
   * @param grid Angular grid to use.
   * @param interaction_variables Interaction variables.
   * @param Tmatrix T-matrix to use.
   * @param nfactors N based resources to use.
   * @param wigner_d Special Wigner D resources to use.
   * @param result Resource in which the result is stored.
   */
  inline ScatteringMatrixTask(
      const ScatteringMatrixGrid &grid,
      const InteractionVariables &interaction_variables,
      const TMatrixResource &Tmatrix, const NBasedResources &nfactors,
      const ScatteringMatrixSpecialWignerDResources &wigner_d,
      ScatteringMatrixResult &result)
      : _grid(grid), _interaction_variables(interaction_variables),
        _Tmatrix(Tmatrix), _nfactors(nfactors), _wigner_d(wigner_d),
        _result(result) {}

  virtual ~ScatteringMatrixTask() {}

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
   * from the given input angles to the given output angles at a particle with
   * its symmetry axis fixed to the @f$z@f$-axis of the reference frame.
   *
   * @param grid_in Internal input grid to sample.
   * @param itheta_in Index of the input zenith angle.
   * @param grid_out Internal output grid to sample.
   * @param itheta_out Index of the output zenith angle.
   * @param cosphi_out Cosine of the output azimuth angle,
   * @f$\cos(\phi{}_s)@f$.
   * @param sinphi_out Sine of the output azimuth angle,
   * @f$\sin(\phi{}_s)@f$.
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<std::complex<float_type>> get_forward_scattering_matrix(
      const uint_fast8_t grid_in, const uint_fast32_t itheta_in,
      const uint_fast8_t grid_out, const uint_fast32_t itheta_out,
      const float_type cosphi_out, const float_type sinphi_out) const {

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
            m * _wigner_d.get_wigner_d_sinx(grid_in, m, itheta_in, nn, nmax);
        const float_type tau_nn =
            _wigner_d.get_dwigner_d(grid_in, m, itheta_in, nn, nmax);

        for (uint_fast32_t n = nmin; n < nmax + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n =
              m * _wigner_d.get_wigner_d_sinx(grid_out, m, itheta_out, n, nmax);
          const float_type tau_n =
              _wigner_d.get_dwigner_d(grid_out, m, itheta_out, n, nmax);

          // get the c factor for these values of n and n'
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

    const uint_fast32_t ntheta_out = _grid._theta_out.size();
    const uint_fast32_t nphi_out = _grid._phi_out.size();

    const std::complex<float_type> icompl(0., 1.);
    const float_type half(0.5);
    for (uint_fast32_t itheta_out = 0; itheta_out < ntheta_out; ++itheta_out) {
      for (uint_fast32_t iphi_out = 0; iphi_out < nphi_out; ++iphi_out) {
        const float_type cos_phi_out = _grid._cos_phi_out[iphi_out];
        const float_type sin_phi_out = _grid._sin_phi_out[iphi_out];

        const Matrix<std::complex<float_type>> S =
            get_forward_scattering_matrix(0, 0, 1, itheta_out, cos_phi_out,
                                          sin_phi_out);

        Matrix<float_type> &Z =
            _result._scattering_matrix[itheta_out * nphi_out + iphi_out];

        Z(0, 0) = (half * (S(0, 0) * conj(S(0, 0)) + S(0, 1) * conj(S(0, 1)) +
                           S(1, 0) * conj(S(1, 0)) + S(1, 1) * conj(S(1, 1))))
                      .real();
        Z(0, 1) = (half * (S(0, 0) * conj(S(0, 0)) - S(0, 1) * conj(S(0, 1)) +
                           S(1, 0) * conj(S(1, 0)) - S(1, 1) * conj(S(1, 1))))
                      .real();
        Z(0, 2) = (-S(0, 0) * conj(S(0, 1)) - S(1, 1) * conj(S(1, 0))).real();
        Z(0, 3) = (icompl * (S(0, 0) * conj(S(0, 1)) - S(1, 1) * conj(S(1, 0))))
                      .real();

        Z(1, 0) = (half * (S(0, 0) * conj(S(0, 0)) + S(0, 1) * conj(S(0, 1)) -
                           S(1, 0) * conj(S(1, 0)) - S(1, 1) * conj(S(1, 1))))
                      .real();
        Z(1, 1) = (half * (S(0, 0) * conj(S(0, 0)) - S(0, 1) * conj(S(0, 1)) -
                           S(1, 0) * conj(S(1, 0)) + S(1, 1) * conj(S(1, 1))))
                      .real();
        Z(1, 2) = (-S(0, 0) * conj(S(0, 1)) + S(1, 1) * conj(S(1, 0))).real();
        Z(1, 3) = (icompl * (S(0, 0) * conj(S(0, 1)) + S(1, 1) * conj(S(1, 0))))
                      .real();

        Z(2, 0) = (-S(0, 0) * conj(S(1, 0)) - S(1, 1) * conj(S(0, 1))).real();
        Z(2, 1) = (-S(0, 0) * conj(S(1, 0)) + S(1, 1) * conj(S(0, 1))).real();
        Z(2, 2) = (S(0, 0) * conj(S(1, 1)) + S(0, 1) * conj(S(1, 0))).real();
        Z(2, 3) =
            (-icompl * (S(0, 0) * conj(S(1, 1)) + S(1, 0) * conj(S(0, 1))))
                .real();

        Z(3, 0) = (icompl * (S(1, 0) * conj(S(0, 0)) + S(1, 1) * conj(S(0, 1))))
                      .real();
        Z(3, 1) = (icompl * (S(1, 0) * conj(S(0, 0)) - S(1, 1) * conj(S(0, 1))))
                      .real();
        Z(3, 2) =
            (-icompl * (S(1, 1) * conj(S(0, 0)) - S(0, 1) * conj(S(1, 0))))
                .real();
        Z(3, 3) = (S(1, 1) * conj(S(0, 0)) - S(0, 1) * conj(S(1, 0))).real();
      }
    }
  }
};

#endif // SCATTERINGMATRIXTASK_HPP
