/**
 * @file AbsorptionCoefficientTask.hpp
 *
 * @brief Task that computes AbsorptionCoefficients.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ABSORPTIONCOEFFICIENTTASK_HPP
#define ABSORPTIONCOEFFICIENTTASK_HPP

#include "Configuration.hpp"
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

  /*! @brief Give access to the averaging task. */
  friend class AbsorptionShapeAveragingTask;

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
 *
 * The grid stores two sets of angles and some useful function values based on
 * these:
 *  - the output zenith angles that correspond to the directions for which we
 *    want to compute absorption coefficients. These are provided by the user
 *    and can be arbitrary.
 *  - the input zenith and azimuth angles that correspond to a grid of incoming
 *    directions over which we integrate. These are computed as the quadrature
 *    points for a Gauss-Legendre quadrature rule. For these we also store the
 *    associated Gauss-Legendre quadrature weights.
 */
class AbsorptionCoefficientGrid : public Resource,
                                  public Task,
                                  public Computable {

  /*! @brief Give access to the computation task. */
  friend class AbsorptionCoefficientTask;

  /*! @brief Give access to special Wigner D resources. */
  friend class AbsorptionSpecialWignerDResources;

private:
  /*! @brief Output zenith angles (in radians). */
  std::vector<float_type> _theta_out;

  /*! @brief Cosines of the output zenith angles. */
  std::vector<float_type> _cos_theta_out;

  /*! @brief Sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out;

  /*! @brief Inverse sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out_inverse;

  /*! @brief Cosines of the input zenith angles. */
  std::vector<float_type> _cos_theta_in;

  /*! @brief Sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in;

  /*! @brief Inverse sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in_inverse;

  /*! @brief Gauss-Legendre weights for cosines of the input zenith angles. */
  std::vector<float_type> _cos_theta_in_weights;

  /*! @brief Cosines of the input azimuth angles. */
  std::vector<float_type> _cos_phi_in;

  /*! @brief Sines of the input azimuth angles. */
  std::vector<float_type> _sin_phi_in;

  /*! @brief Gauss-Legendre weights for the input azimuth angles. */
  std::vector<float_type> _phi_in_weights;

public:
  /**
   * @brief Constructor.
   *
   * @param ntheta Number of zenith angles.
   * @param theta Grid of zenith angles (in radians, of size ntheta or more).
   * @param ngauss Number of Gauss-Legendre quadrature points for directional
   * averaging.
   */
  inline AbsorptionCoefficientGrid(const uint_fast32_t ntheta,
                                   const float_type *theta,
                                   const uint_fast32_t ngauss)
      : _theta_out(ntheta), _cos_theta_out(ntheta), _sin_theta_out(ntheta),
        _sin_theta_out_inverse(ntheta), _cos_theta_in(ngauss),
        _sin_theta_in(ngauss), _sin_theta_in_inverse(ngauss),
        _cos_theta_in_weights(ngauss), _cos_phi_in(ngauss), _sin_phi_in(ngauss),
        _phi_in_weights(ngauss) {

    for (uint_fast32_t i = 0; i < ntheta; ++i) {
      _theta_out[i] = theta[i];
    }
  }

  virtual ~AbsorptionCoefficientGrid() {}

  /**
   * @brief Get the size in memory of a hypothetical AbsorptionCoefficientGrid
   * object with the given parameters.
   *
   * @param ntheta Number of zenith angles.
   * @param ngauss Number of Gauss-Legendre quadrature points for directional
   * averaging.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t ntheta,
                                       const uint_fast32_t ngauss) {
    size_t size = sizeof(AbsorptionCoefficientGrid);
    size += 4 * ntheta * sizeof(float_type);
    size += 6 * ngauss * sizeof(float_type);
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
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 0.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 0; }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _cos_theta_out.size();
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      ctm_assert(_theta_out[itheta] >= 0.);
      ctm_assert(_theta_out[itheta] <= M_PI);
      _cos_theta_out[itheta] = cos(_theta_out[itheta]);
      // note that we only use the positive root since the sin(theta) values
      // are guaranteed to be in the range [0, 1]
      _sin_theta_out[itheta] =
          sqrt((1. - _cos_theta_out[itheta]) * (1. + _cos_theta_out[itheta]));
      if (_sin_theta_out[itheta] != 0.) {
        _sin_theta_out_inverse[itheta] = 1. / _sin_theta_out[itheta];
      } else {
        // this value should be ignored, we just set it to a safe and
        // recognisable value
        _sin_theta_out_inverse[itheta] = 9000.;
      }
    }

    const uint_fast32_t ngauss = _cos_theta_in.size();
    SpecialFunctions::get_gauss_legendre_points_and_weights<float_type>(
        ngauss, _cos_theta_in, _cos_theta_in_weights);
    // note that we temporarily store the phi angles in cos_phi
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        ngauss, 0., 2. * M_PI, _cos_phi_in, _phi_in_weights);

    for (uint_fast32_t igauss = 0; igauss < ngauss; ++igauss) {
      // again, we know that the sin(theta) values are in the range [0, 1]
      _sin_theta_in[igauss] =
          sqrt((1. - _cos_theta_in[igauss]) * (1. + _cos_theta_in[igauss]));
      if (_sin_theta_in[igauss] != 0.) {
        _sin_theta_in_inverse[igauss] = 1. / _sin_theta_in[igauss];
      } else {
        // again: safe and recognisable value
        _sin_theta_in_inverse[igauss] = 9000.;
      }

      // convert the phi angles to cos_phi
      const float_type phi_out = _cos_phi_in[igauss];
      _cos_phi_in[igauss] = cos(phi_out);
      // note that in this case we cannot use a sqrt(1-cos^2) rule, since we
      // can have sin(phi)<0.
      _sin_phi_in[igauss] = sin(phi_out);
    }
  }

  /**
   * @brief Get the number of angles in the grid.
   *
   * @return Number of angles.
   */
  inline uint_fast32_t get_number_of_angles() const {
    return _theta_out.size();
  }
};

/**
 * @brief Precomputed special Wigner D functions that depend on a specific value
 * of @f$n_{max}@f$ and @f$n_{GL}@f$.
 */
class AbsorptionSpecialWignerDResources : public Resource,
                                          public Task,
                                          public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Absorption coefficient grid to use. */
  const AbsorptionCoefficientGrid &_grid;

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
  inline AbsorptionSpecialWignerDResources(
      const uint_fast32_t nmax, const AbsorptionCoefficientGrid &grid)
      : _nmax(nmax), _grid(grid) {

    const uint_fast32_t number_of_elements = (2 + nmax + 1) * nmax / 2;

    _wigner_d_sinx[0].reserve(number_of_elements);
    _wigner_d_sinx[1].reserve(number_of_elements);
    _dwigner_d[0].reserve(number_of_elements);
    _dwigner_d[1].reserve(number_of_elements);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        _wigner_d_sinx[0].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
        _dwigner_d[0].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
        _wigner_d_sinx[1].push_back(
            Matrix<float_type>(grid._cos_theta_in.size(), n));
        _dwigner_d[1].push_back(
            Matrix<float_type>(grid._cos_theta_in.size(), n));
      }
    }
    ctm_assert(_wigner_d_sinx[0].size() == number_of_elements);
  }

  virtual ~AbsorptionSpecialWignerDResources() {}

  /**
   * @brief Get the size in memory of a hypothetical SpecialWignerDResources
   * object with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ntheta Number of zenith angles.
   * @return Size in bytes of the hypothetical object.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ntheta) {
    size_t size = sizeof(AbsorptionSpecialWignerDResources);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        size += 2 * ntheta * n * sizeof(float_type);
        size += 2 * ntheta * n * sizeof(float_type);
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
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 1; }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    const uint_fast32_t ntheta_in = _grid._cos_theta_out.size();
    const uint_fast32_t ntheta_out = _grid._cos_theta_in.size();
    for (uint_fast32_t n = 1; n < _nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        for (uint_fast32_t itheta_in = 0; itheta_in < ntheta_in; ++itheta_in) {
          const uint_fast32_t index = (2 + n) * (n - 1) / 2 + m;
          ctm_assert(index < _wigner_d_sinx[0].size());
          ctm_assert_message(
              _wigner_d_sinx[0][index].get_number_of_columns() == n,
              "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
              ", index: %" PRIuFAST32,
              _wigner_d_sinx[0][index].get_number_of_columns(), n, m, index);
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_out[itheta_in], _grid._sin_theta_out[itheta_in],
              _grid._sin_theta_out_inverse[itheta_in], n, m,
              &_wigner_d_sinx[0][index].get_row(itheta_in)[0],
              &_dwigner_d[0][index].get_row(itheta_in)[0]);
        }
        for (uint_fast32_t itheta_out = 0; itheta_out < ntheta_out;
             ++itheta_out) {
          const uint_fast32_t index = (2 + n) * (n - 1) / 2 + m;
          ctm_assert(index < _wigner_d_sinx[1].size());
          ctm_assert_message(
              _wigner_d_sinx[1][index].get_number_of_columns() == n,
              "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
              ", index: %" PRIuFAST32,
              _wigner_d_sinx[1][index].get_number_of_columns(), n, m, index);
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_in[itheta_out], _grid._sin_theta_in[itheta_out],
              _grid._sin_theta_in_inverse[itheta_out], n, m,
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
 * @brief Task that computes AbsorptionCoefficients.
 */
class AbsorptionCoefficientTask : public Task {
private:
  /*! @brief Zenith angle grid to use. */
  const AbsorptionCoefficientGrid &_grid;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief T-matrix to use (read only). */
  const TMatrixResource &_Tmatrix;

  /*! @brief N based resources to use (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Special Wigner D resources to use (read only). */
  const AbsorptionSpecialWignerDResources &_wigner_d;

  /*! @brief Resource in which the result is stored. */
  AbsorptionCoefficientResult &_result;

  /*! @brief Subtract the directionally averaged scattering cross sections? */
  const bool _account_for_scattering;

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
   * @param account_for_scattering Subtract the directionally averaged
   * scattering cross sections?
   */
  inline AbsorptionCoefficientTask(
      const AbsorptionCoefficientGrid &grid,
      const InteractionVariables &interaction_variables,
      const TMatrixResource &Tmatrix, const NBasedResources &nfactors,
      const AbsorptionSpecialWignerDResources &wigner_d,
      AbsorptionCoefficientResult &result,
      const bool account_for_scattering = false)
      : _grid(grid), _interaction_variables(interaction_variables),
        _Tmatrix(Tmatrix), _nfactors(nfactors), _wigner_d(wigner_d),
        _result(result), _account_for_scattering(account_for_scattering) {}

  virtual ~AbsorptionCoefficientTask() {}

  /**
   * @brief Get the size in memory of a hypothetical AbsorptionCoefficientTask
   * object with the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(AbsorptionCoefficientTask);
    return size;
  }

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
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 5.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 5; }

  /**
   * @brief Get the forward scattering matrix @f$S@f$ for a scattering event
   * from the given input angles to the given output angles at a particle with
   * its symmetry axis fixed to the @f$z@f$-axis of the reference frame.
   *
   * @param grid_out Internal output grid to sample.
   * @param itheta_out Index of the output zenith angle.
   * @param grid_in Internal input grid to sample.
   * @param itheta_in Index of the input zenith angle.
   * @param cosphi_in Cosine of the input azimuth angle, @f$\cos(\phi{}_i)@f$.
   * @param sinphi_in Sine of the input azimuth angle, @f$\sin(\phi{}_i)@f$.
   * @return Scattering matrix for this scattering event.
   */
  inline Matrix<std::complex<float_type>> get_forward_scattering_matrix(
      const uint_fast8_t grid_in, const uint_fast32_t itheta_in,
      const uint_fast8_t grid_out, const uint_fast32_t itheta_out,
      const float_type cosphi_in, const float_type sinphi_in) const {

    const uint_fast32_t nmax = _Tmatrix.get_nmax();

    // now compute the matrix S^P
    // we precompute e^{i(phi_out-phi_in)}
    const std::complex<float_type> expiphi_p_out_m_in(cosphi_in, -sinphi_in);
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

    const uint_fast32_t ntheta = _grid._theta_out.size();
    const uint_fast32_t ngauss = _grid._cos_theta_in.size();

    for (uint_fast32_t itheta_out = 0; itheta_out < ntheta; ++itheta_out) {

      Matrix<std::complex<float_type>> S =
          get_forward_scattering_matrix(0, itheta_out, 0, itheta_out, 1., 0.);

      const float_type prefactor =
          2. * M_PI / _interaction_variables.get_wavenumber();
      _result._Qabs[itheta_out] = prefactor * (S(0, 0) + S(1, 1)).imag();
      _result._Qabspol[itheta_out] = prefactor * (S(0, 0) - S(1, 1)).imag();

      if (_account_for_scattering) {
        const float_type half(0.5);
        for (uint_fast32_t itheta_in = 0; itheta_in < ngauss; ++itheta_in) {
          for (uint_fast32_t iphi_in = 0; iphi_in < ngauss; ++iphi_in) {
            const float_type cos_phi_in = _grid._cos_phi_in[iphi_in];
            const float_type sin_phi_in = _grid._sin_phi_in[iphi_in];

            Matrix<std::complex<float_type>> Stp =
                get_forward_scattering_matrix(1, itheta_in, 0, itheta_out,
                                              cos_phi_in, sin_phi_in);

            const float_type weight = _grid._cos_theta_in_weights[itheta_in] *
                                      _grid._phi_in_weights[iphi_in];
            const float_type Z00 =
                (half *
                 (Stp(0, 0) * conj(Stp(0, 0)) + Stp(0, 1) * conj(Stp(0, 1)) +
                  Stp(1, 0) * conj(Stp(1, 0)) + Stp(1, 1) * conj(Stp(1, 1))))
                    .real();
            _result._Qabs[itheta_out] -= Z00 * weight;

            const float_type Z10 =
                (half *
                 (Stp(0, 0) * conj(Stp(0, 0)) + Stp(0, 1) * conj(Stp(0, 1)) -
                  Stp(1, 0) * conj(Stp(1, 0)) - Stp(1, 1) * conj(Stp(1, 1))))
                    .real();
            _result._Qabspol[itheta_out] -= Z10 * weight;
          }
        }
      }

      // normalise the coefficients
      const float_type a = _interaction_variables.get_equal_volume_radius();
      const float_type norm = M_PI * a * a;
      _result._Qabs[itheta_out] /= norm;
      _result._Qabspol[itheta_out] /= norm;
    }
  }
};

#endif // ABSORPTIONCOEFFICIENTTASK_HPP
