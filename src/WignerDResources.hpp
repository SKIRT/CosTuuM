/**
 * @file WignerDResources.hpp
 *
 * @brief Precomputed Wigner D functions that depend on a specific value of
 * @f$n_{max}@f$ and @f$n_{GL}@f$.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef WIGNERDRESOURCES_HPP
#define WIGNERDRESOURCES_HPP

#include "Configuration.hpp"
#include "ConvergedSizeResources.hpp"
#include "GaussBasedResources.hpp"
#include "Matrix.hpp"
#include "QuickSchedWrapper.hpp"
#include "SpecialFunctions.hpp"

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Precomputed Wigner D functions that depend on a specific value of
 * @f$n_{max}@f$ and @f$n_{GL}@f$.
 */
class WignerDm0Resources : public Resource, public Task, public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Wigner D functions. */
  Matrix<float_type> _wigner_d;

  /*! @brief Derivatives of the Wigner D functions. */
  Matrix<float_type> _dwigner_d;

  /*! @brief Gauss-Legendre quadrature points. */
  const GaussBasedResources &_quadrature_points;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param quadrature_points Gauss-Legendre quadrature points (to be computed
   * before this task is executed!).
   */
  inline WignerDm0Resources(const uint_fast32_t nmax,
                            const uint_fast32_t ngauss,
                            const GaussBasedResources &quadrature_points)
      : _nmax(nmax), _ngauss(ngauss), _wigner_d(2 * ngauss, nmax),
        _dwigner_d(2 * ngauss, nmax), _quadrature_points(quadrature_points) {}

  virtual ~WignerDm0Resources() {}

  /**
   * @brief Get the size in memory of a hypothetical WignerDResources object
   * with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ngauss) {
    size_t size = sizeof(WignerDm0Resources);
    // wigner_d
    size += 2 * ngauss * nmax * sizeof(float_type);
    // dwigner_d
    size += 2 * ngauss * nmax * sizeof(float_type);
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
    quicksched.link_task_and_resource(*this, _quadrature_points, false);
  }

  /**
   * @brief Compute the factors.
   */
  virtual void execute() {
    for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {
      const uint_fast32_t i1 = _ngauss + ig;
      const uint_fast32_t i2 = _ngauss - ig + 1;
      std::vector<float_type> dv1(_nmax), dv2(_nmax);
      SpecialFunctions::wigner_dn_0m(_quadrature_points.get_costheta(i1 - 1),
                                     _nmax, 0, &dv1[0], &dv2[0]);
      int_fast8_t sign = 1;
      for (uint_fast32_t n = 0; n < _nmax; ++n) {
        sign = -sign;
        _wigner_d(i1 - 1, n) = dv1[n];
        _wigner_d(i2 - 1, n) = sign * dv1[n];
        _dwigner_d(i1 - 1, n) = dv2[n];
        _dwigner_d(i2 - 1, n) = -sign * dv2[n];

        ctm_assert_not_nan(_wigner_d(i1 - 1, n));
        ctm_assert_not_nan(_wigner_d(i2 - 1, n));
        ctm_assert_not_nan(_dwigner_d(i1 - 1, n));
        ctm_assert_not_nan(_dwigner_d(i2 - 1, n));
      }
    }
    make_available();
  }

  /**
   * @brief Get the Wigner D function for the given Gauss-Legendre quadrature
   * point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function value.
   */
  inline float_type get_wigner_d(const uint_fast32_t ig,
                                 const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();
    return _wigner_d(ig, n - 1);
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given Gauss-
   * Legendre quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast32_t ig,
                                  const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();
    return _dwigner_d(ig, n - 1);
  }
};

/**
 * @brief Precomputed Wigner D functions that depend on a specific value of
 * @f$n_{max}@f$ and @f$n_{GL}@f$ that is read from a TMatrixResource.
 */
class WignerDmn0Resources : public Resource, public Task, public Computable {
private:
  /*! @brief @f$m@f$ value for which the factors are computed. */
  const uint_fast32_t _m;

  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Maximum number of Gauss-Legendre quadrature points,
   *  @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Wigner D functions. */
  Matrix<float_type> _wigner_d;

  /*! @brief Derivatives of the Wigner D functions. */
  Matrix<float_type> _dwigner_d;

  /*! @brief Gauss-Legendre quadrature points. */
  const GaussBasedResources &_quadrature_points;

  /*! @brief Actual order and number of quadrature points. */
  const ConvergedSizeResources &_converged_size;

public:
  /**
   * @brief Constructor.
   *
   * @param m @f$m@f$ value for which the factors are computed.
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Maximum number of Gauss-Legendre quadrature points,
   * @f$n_{GL}@f$.
   * @param quadrature_points Gauss-Legendre quadrature points (to be computed
   * before this task is executed!).
   * @param converged_size Actual order and number of quadrature points.
   */
  inline WignerDmn0Resources(const uint_fast32_t m, const uint_fast32_t nmax,
                             const uint_fast32_t ngauss,
                             const GaussBasedResources &quadrature_points,
                             const ConvergedSizeResources &converged_size)
      : _m(m), _nmax(nmax), _ngauss(ngauss), _wigner_d(2 * ngauss, nmax),
        _dwigner_d(2 * ngauss, nmax), _quadrature_points(quadrature_points),
        _converged_size(converged_size) {}

  virtual ~WignerDmn0Resources() {}

  /**
   * @brief Get the size in memory of a hypothetical WignerDResources object
   * with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ngauss) {
    size_t size = sizeof(WignerDm0Resources);
    // wigner_d
    size += 2 * ngauss * nmax * sizeof(float_type);
    // dwigner_d
    size += 2 * ngauss * nmax * sizeof(float_type);
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
    quicksched.link_task_and_resource(*this, _quadrature_points, false);
    quicksched.link_task_and_resource(*this, _converged_size, false);
  }

  /**
   * @brief Compute the factors.
   */
  virtual void execute() {

    const uint_fast32_t nmax = _converged_size.get_nmax();
    const uint_fast32_t ngauss = _converged_size.get_ngauss();

    ctm_assert(nmax <= _nmax);
    ctm_assert(ngauss <= _ngauss);

    for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
      const uint_fast32_t i1 = ngauss + ig;
      const uint_fast32_t i2 = ngauss - ig + 1;
      std::vector<float_type> dv1(nmax), dv2(nmax);
      SpecialFunctions::wigner_dn_0m(_quadrature_points.get_costheta(i1 - 1),
                                     nmax, _m, &dv1[0], &dv2[0]);
      int_fast8_t sign = 1;
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        sign = -sign;
        _wigner_d(i1 - 1, n) = dv1[n];
        _wigner_d(i2 - 1, n) = sign * dv1[n];
        _dwigner_d(i1 - 1, n) = dv2[n];
        _dwigner_d(i2 - 1, n) = -sign * dv2[n];

        ctm_assert_not_nan(_wigner_d(i1 - 1, n));
        ctm_assert_not_nan(_wigner_d(i2 - 1, n));
        ctm_assert_not_nan(_dwigner_d(i1 - 1, n));
        ctm_assert_not_nan(_dwigner_d(i2 - 1, n));
      }
    }
    make_available();
  }

  /**
   * @brief Get the Wigner D function for the given Gauss-Legendre quadrature
   * point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function value.
   */
  inline float_type get_wigner_d(const uint_fast32_t ig,
                                 const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();
    return _wigner_d(ig, n - 1);
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given Gauss-
   * Legendre quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast32_t ig,
                                  const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();
    return _dwigner_d(ig, n - 1);
  }
};

#endif // WIGNERDRESOURCES_HPP
