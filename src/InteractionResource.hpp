/**
 * @file InteractionResource.hpp
 *
 * @brief Precomputed factors that depend on precomputed Gauss-Legendre factors,
 * precomputed, particle specific geometric factors, and a specific choice of
 * particle refractive index and interaction wavelength.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef INTERACTIONRESOURCE_HPP
#define INTERACTIONRESOURCE_HPP

#include "Configuration.hpp"
#include "ConvergedSizeResources.hpp"
#include "GaussBasedResources.hpp"
#include "Matrix.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"

#include <complex>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Scalar variables involved in a particle - photon interaction.
 */
class InteractionVariables {
private:
  /*! @brief Equal volume sphere radius, @f$R_V@f$ (in m). */
  const float_type _R_V;

  /*! @brief Wavenumber, @f$k = \frac{2\pi{}}{\lambda{}}@f$. */
  const float_type _k;

  /*! @brief Wavenumber squared, @f$k^2@f$. */
  const float_type _k2;

  /*! @brief Wavenumber times refractive index, @f$m_rk@f$. */
  const std::complex<float_type> _kmr;

  /*! @brief Wavenumber squared times refractive index, @f$m_rk^2@f$. */
  const std::complex<float_type> _k2mr;

public:
  /**
   * @brief Constructor.
   *
   * @param R_V Equal volume sphere radius, @f$R_V@f$ (in m).
   * @param wavelength Wavelength of incident radiation, @f$\lambda{}@f$ (in m).
   * @param refractive_index Refractive index of the material, @f$m_r@f$.
   */
  inline InteractionVariables(const float_type R_V, const float_type wavelength,
                              const std::complex<float_type> refractive_index)
      : _R_V(R_V), _k(2. * M_PI / wavelength), _k2(_k * _k),
        _kmr(refractive_index * _k), _k2mr(refractive_index * _k2) {}

  /**
   * @brief Get the equal volume radius, @f$R_V@f$.
   *
   * @return Equal volume radius, @f$R_V@f$.
   */
  inline float_type get_equal_volume_radius() const { return _R_V; }

  /**
   * @brief Get the wavenumber, @f$k@f$.
   *
   * @return Wavenumber, @f$k@f$.
   */
  inline float_type get_wavenumber() const { return _k; }

  /**
   * @brief Get the wavenumber squared, @f$k^2@f$.
   *
   * @return Wavenumber squared, @f$k^2@f$.
   */
  inline float_type get_wavenumber_squared() const { return _k2; }

  /**
   * @brief Get the wavenumber times the refractive index of the material,
   * @f$m_rk@f$.
   *
   * @return Wavenumber times refractive index, @f$m_rk@f$.
   */
  inline std::complex<float_type> get_material_wavenumber() const {
    return _kmr;
  }

  /**
   * @brief Get the wavenumber squared times the refractive index of the
   * material, @f$m_rk^2@f$.
   *
   * @return Wavenumber squared times refractive index, @f$m_rk^2@f$.
   */
  inline std::complex<float_type>
  get_material_wavenumber_times_wavenumber() const {
    return _k2mr;
  }
};

/**
 * @brief Precomputed factors that depend on precomputed Gauss-Legendre factors,
 * precomputed, particle specific geometric factors, and a specific choice of
 * particle refractive index and interaction wavelength.
 */
class InteractionResource : public Resource {

  /*! @brief Grant access to computation tasks. */
  friend class InteractionTask;

private:
  /*! @brief Precomputed factors @f$kr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<float_type> _kr;

  /*! @brief Precomputed factors @f$\frac{1}{kr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _krinv;

  /*! @brief Precomputed factors @f$km_rr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<std::complex<float_type>> _krmr;

  /*! @brief Precomputed factors @f$\frac{1}{km_rr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<std::complex<float_type>> _krmrinv;

  /*! @brief Bessel functions @f$j_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<float_type> _jkr;

  /*! @brief Bessel functions @f$y_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<float_type> _ykr;

  /*! @brief Bessel function derivatives @f$\frac{[krj_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _djkr;

  /*! @brief Bessel function derivatives @f$\frac{[kry_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _dykr;

  /*! @brief Bessel functions @f$j_n(km_rr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<float_type>> _jkrmr;

  /*! @brief Bessel function derivatives
   *  @f$\frac{[km_rrj(km_rr)]'}{km_rr}@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<float_type>> _djkrmr;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order of spherical basis, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline InteractionResource(const uint_fast32_t nmax,
                             const uint_fast32_t ngauss)
      : _kr(2 * ngauss, float_type(0.)), _krinv(2 * ngauss, float_type(0.)),
        _krmr(2 * ngauss, float_type(0.)), _krmrinv(2 * ngauss, float_type(0.)),
        _jkr(2 * ngauss, nmax), _ykr(2 * ngauss, nmax), _djkr(2 * ngauss, nmax),
        _dykr(2 * ngauss, nmax), _jkrmr(2 * ngauss, nmax),
        _djkrmr(2 * ngauss, nmax) {}

  virtual ~InteractionResource() {}

  /**
   * @brief Get the size in memory of a hypothetical InteractionResource object
   * with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ngauss) {
    size_t size = sizeof(InteractionResource);
    // kr
    size += 2 * ngauss * sizeof(float_type);
    // krinv
    size += 2 * ngauss * sizeof(float_type);
    // krmr
    size += 4 * ngauss * sizeof(float_type);
    // krmrinv
    size += 4 * ngauss * sizeof(float_type);
    // jkr
    size += 2 * ngauss * nmax * sizeof(float_type);
    // ykr
    size += 2 * ngauss * nmax * sizeof(float_type);
    // djkr
    size += 2 * ngauss * nmax * sizeof(float_type);
    // dykr
    size += 2 * ngauss * nmax * sizeof(float_type);
    // jkrmr
    size += 4 * ngauss * nmax * sizeof(float_type);
    // djkrmr
    size += 4 * ngauss * nmax * sizeof(float_type);
    return size;
  }

  /**
   * @brief Get the computational cost of this task.
   *
   * @return Computational cost.
   */
  virtual int_fast32_t get_cost() const { return 157675 * _kr.size() - 755642; }

  /**
   * @brief Get the wavenumber multiplied with the radius for the given Gauss-
   * Legendre quadrature point.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @return Corresponding value of @f$kr@f$.
   */
  inline float_type get_kr(const uint_fast32_t ig) const {
    ctm_assert(ig < _kr.size());
    return _kr[ig];
  }

  /**
   * @brief Get @f$\frac{1}{kr}@f$ for the given Gauss-Legendre quadrature
   * point.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @return Corresponding value of @f$\frac{1}{kr}@f$.
   */
  inline float_type get_krinv(const uint_fast32_t ig) const {
    ctm_assert(ig < _krinv.size());
    return _krinv[ig];
  }

  /**
   * @brief Get the wavenumber multiplied with the radius and refractive index
   * for the given Gauss-Legendre quadrature point.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @return Corresponding value of @f$km_rr@f$.
   */
  inline std::complex<float_type> get_krmr(const uint_fast32_t ig) const {
    ctm_assert(ig < _krmr.size());
    return _krmr[ig];
  }

  /**
   * @brief Get @f$\frac{1}{km_rr}@f$ for the given Gauss-Legendre quadrature
   * point.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @return Corresponding value of @f$\frac{1}{km_rr}@f$.
   */
  inline std::complex<float_type> get_krmrinv(const uint_fast32_t ig) const {
    ctm_assert(ig < _krmrinv.size());
    return _krmrinv[ig];
  }

  /**
   * @brief Get @f$j_n(kr)@f$ for the given Gauss-Legendre quadrature point and
   * order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$j_n(kr)@f$.
   */
  inline float_type get_jkr(const uint_fast32_t ig,
                            const uint_fast32_t n) const {
    ctm_assert(ig < _jkr.get_number_of_rows());
    ctm_assert(n - 1 < _jkr.get_number_of_columns());
    return _jkr(ig, n - 1);
  }

  /**
   * @brief Get @f$y_n(kr)@f$ for the given Gauss-Legendre quadrature point and
   * order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$y_n(kr)@f$.
   */
  inline float_type get_ykr(const uint_fast32_t ig,
                            const uint_fast32_t n) const {
    ctm_assert(ig < _ykr.get_number_of_rows());
    ctm_assert(n - 1 < _ykr.get_number_of_columns());
    return _ykr(ig, n - 1);
  }

  /**
   * @brief Get @f$\frac{[krj_n(kr)]'}{kr}@f$ for the given Gauss-Legendre
   * quadrature point and order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$\frac{[krj_n(kr)]'}{kr}@f$
   */
  inline float_type get_djkr(const uint_fast32_t ig,
                             const uint_fast32_t n) const {
    ctm_assert(ig < _djkr.get_number_of_rows());
    ctm_assert(n - 1 < _djkr.get_number_of_columns());
    return _djkr(ig, n - 1);
  }

  /**
   * @brief Get @f$\frac{[kry_n(kr)]'}{kr}@f$ for the given Gauss-Legendre
   * quadrature point and order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$\frac{[kry_n(kr)]'}{kr}@f$
   */
  inline float_type get_dykr(const uint_fast32_t ig,
                             const uint_fast32_t n) const {
    ctm_assert(ig < _dykr.get_number_of_rows());
    ctm_assert(n - 1 < _dykr.get_number_of_columns());
    return _dykr(ig, n - 1);
  }

  /**
   * @brief Get @f$j_n(km_rr)@f$ for the given Gauss-Legendre quadrature point
   * and order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$j_n(km_rr)@f$.
   */
  inline std::complex<float_type> get_jkrmr(const uint_fast32_t ig,
                                            const uint_fast32_t n) const {
    ctm_assert(ig < _jkrmr.get_number_of_rows());
    ctm_assert(n - 1 < _jkrmr.get_number_of_columns());
    return _jkrmr(ig, n - 1);
  }

  /**
   * @brief Get @f$\frac{[km_rrj_n(km_rr)]'}{km_rr}@f$ for the given Gauss-
   * Legendre quadrature point and order.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return @f$\frac{[km_rrj_n(km_rr)]'}{km_rr}@f$
   */
  inline std::complex<float_type> get_djkrmr(const uint_fast32_t ig,
                                             const uint_fast32_t n) const {
    ctm_assert(ig < _djkrmr.get_number_of_rows());
    ctm_assert(n - 1 < _djkrmr.get_number_of_columns());
    return _djkrmr(ig, n - 1);
  }
};

/**
 * @brief Task that computes interaction resources for a specific order and
 * number of quadrature points.
 */
class InteractionTask : public Task {
private:
  /*! @brief Maximum order of the spherical basis functions, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Geometrical factors for this particle (read only). */
  const ParticleGeometryResource &_geometry;

  /*! @brief Resource keeping track of the convergence of the calculation (read
   *  only). */
  const ConvergedSizeResources &_converged_size;

  /*! @brief Interaction variables for this interaction (read only). */
  const InteractionVariables &_interaction_variables;

  /*! @brief InteractionResource object to operate on. */
  InteractionResource &_interaction;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order of the spherical basis functions, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param geometry Geometrical factors for this particle.
   * @param converged_size Resource keeping track of the convergence of the
   * calculation.
   * @param interaction_variables Interaction variables for this interaction.
   * @param interaction InteractionResource object to operate on.
   */
  inline InteractionTask(const uint_fast32_t nmax, const uint_fast32_t ngauss,
                         const ParticleGeometryResource &geometry,
                         const ConvergedSizeResources &converged_size,
                         const InteractionVariables &interaction_variables,
                         InteractionResource &interaction)
      : _nmax(nmax), _ngauss(ngauss), _geometry(geometry),
        _converged_size(converged_size),
        _interaction_variables(interaction_variables),
        _interaction(interaction) {}

  virtual ~InteractionTask() {}

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _interaction, true);

    // read access
    quicksched.link_task_and_resource(*this, _geometry, false);
    quicksched.link_task_and_resource(*this, _converged_size, false);
  }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    if (_converged_size.is_converged()) {
      // skip this task if the T-matrix was already converged
      return;
    }

    for (uint_fast32_t ig = 0; ig < 2 * _ngauss; ++ig) {
      _interaction._kr[ig] = _interaction_variables.get_wavenumber() *
                             _geometry.get_r(ig) *
                             _interaction_variables.get_equal_volume_radius();
      _interaction._krmr[ig] =
          _interaction_variables.get_material_wavenumber() *
          _geometry.get_r(ig) *
          _interaction_variables.get_equal_volume_radius();
      _interaction._krinv[ig] = float_type(1.) / _interaction._kr[ig];
      _interaction._krmrinv[ig] =
          std::complex<float_type>(1.) / _interaction._krmr[ig];
    }

    ctm_assert_no_nans(_interaction._kr, 2 * _ngauss);
    ctm_assert_no_nans(_interaction._krmr, 2 * _ngauss);
    ctm_assert_no_nans(_interaction._krinv, 2 * _ngauss);
    ctm_assert_no_nans(_interaction._krmrinv, 2 * _ngauss);

    for (uint_fast32_t ig = 0; ig < 2 * _ngauss; ++ig) {
      SpecialFunctions::spherical_j_jdj_array(_nmax, _interaction._kr[ig],
                                              _interaction._jkr.get_row(ig),
                                              _interaction._djkr.get_row(ig));
      SpecialFunctions::spherical_y_ydy_array(_nmax, _interaction._kr[ig],
                                              _interaction._ykr.get_row(ig),
                                              _interaction._dykr.get_row(ig));
      SpecialFunctions::spherical_j_jdj_array(_nmax, _interaction._krmr[ig],
                                              _interaction._jkrmr.get_row(ig),
                                              _interaction._djkrmr.get_row(ig));

      ctm_assert_message_no_nans(_interaction._jkr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
      ctm_assert_message_no_nans(_interaction._djkr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
      ctm_assert_message_no_nans(_interaction._ykr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
      ctm_assert_message_no_nans(_interaction._dykr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
      ctm_assert_message_no_nans(_interaction._jkrmr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
      ctm_assert_message_no_nans(_interaction._djkrmr.get_row(ig), _nmax,
                                 "_nmax: %" PRIuFAST32, _nmax);
    }
  }
};

#endif // INTERACTIONRESOURCE_HPP
