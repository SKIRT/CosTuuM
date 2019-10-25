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
 * @brief Precomputed factors that depend on precomputed Gauss-Legendre factors,
 * precomputed, particle specific geometric factors, and a specific choice of
 * particle refractive index and interaction wavelength.
 */
class InteractionResource : public Resource, public Task {
private:
  /*! @brief Maximum order of the spherical basis functions, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Wavenumber, @f$k = \frac{2\pi{}}{\lambda{}}@f$. */
  const float_type _k;

  /*! @brief Wavenumber squared, @f$k^2@f$. */
  const float_type _k2;

  /*! @brief Wavenumber times refractive index, @f$m_rk@f$. */
  const std::complex<float_type> _kmr;

  /*! @brief Wavenumber squared times refractive index, @f$m_rk^2@f$. */
  const std::complex<float_type> _k2mr;

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

  /*! @brief Geometrical factors for this particle (read only). */
  const ParticleGeometryResource &_geometry;

public:
  /**
   * @brief Constructor.
   *
   * @param wavelength Wavelength of incident radiation, @f$\lambda{}@f$.
   * @param refractive_index Refractive index of the material, @f$m_r@f$.
   * @param nmax Maximum order of spherical basis, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param geometry Geometrical factors for this particle (read only).
   */
  inline InteractionResource(const float_type wavelength,
                             const std::complex<float_type> refractive_index,
                             const uint_fast32_t nmax,
                             const uint_fast32_t ngauss,
                             const ParticleGeometryResource &geometry)
      : _nmax(nmax), _ngauss(ngauss), _k(2. * M_PI / wavelength), _k2(_k * _k),
        _kmr(refractive_index * _k), _k2mr(refractive_index * _k2),
        _kr(2 * ngauss, float_type(0.)), _krinv(2 * ngauss, float_type(0.)),
        _krmr(2 * ngauss, float_type(0.)), _krmrinv(2 * ngauss, float_type(0.)),
        _jkr(2 * ngauss, nmax), _ykr(2 * ngauss, nmax), _djkr(2 * ngauss, nmax),
        _dykr(2 * ngauss, nmax), _jkrmr(2 * ngauss, nmax),
        _djkrmr(2 * ngauss, nmax), _geometry(geometry) {}

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
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);

    // read access
    quicksched.link_task_and_resource(*this, _geometry, false);
  }

  /**
   * @brief Compute the factors.
   */
  virtual void execute() {
    for (uint_fast32_t ig = 0; ig < 2 * _ngauss; ++ig) {
      _kr[ig] = _k * _geometry.get_r(ig);
      _krmr[ig] = _kmr * _geometry.get_r(ig);
      _krinv[ig] = 1. / _kr[ig];
      _krmrinv[ig] = 1. / _krmr[ig];
    }
    for (uint_fast32_t ig = 0; ig < 2 * _ngauss; ++ig) {
      SpecialFunctions::spherical_j_jdj_array(_nmax, _kr[ig], _jkr.get_row(ig),
                                              _djkr.get_row(ig));
      SpecialFunctions::spherical_y_ydy_array(_nmax, _kr[ig], _ykr.get_row(ig),
                                              _dykr.get_row(ig));
      SpecialFunctions::spherical_j_jdj_array(
          _nmax, _krmr[ig], _jkrmr.get_row(ig), _djkrmr.get_row(ig));
    }
  }

  /**
   * @brief Get the wavenumber.
   *
   * @return Wavenumber.
   */
  inline float_type get_k() const { return _k; }

  /**
   * @brief Get the wavenumber squared.
   *
   * @return Wavenumber squared.
   */
  inline float_type get_k2() const { return _k2; }

  /**
   * @brief Get the wavenumber squared times the refractive index.
   *
   * @return Wavenumber squared times the refractive index.
   */
  inline std::complex<float_type> get_k2mr() const { return _k2mr; }

  /**
   * @brief Get the wavenumber multiplied with the radius for the given Gauss-
   * Legendre quadrature point.
   *
   * @param ig Index of a Gauss-Legendre quadrature point.
   * @return Corresponding value of @f$kr@f$.
   */
  inline float_type get_kr(const uint_fast32_t ig) const {
    ctm_assert(ig < 2 * _ngauss);
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
    ctm_assert(ig < 2 * _ngauss);
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
    ctm_assert(ig < 2 * _ngauss);
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
    ctm_assert(ig < 2 * _ngauss);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
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
    ctm_assert(ig < 2 * _ngauss);
    ctm_assert(n - 1 < _nmax);
    return _djkrmr(ig, n - 1);
  }
};

#endif // INTERACTIONRESOURCE_HPP
