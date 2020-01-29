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
class WignerDResources : public Resource, public Task, public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Wigner D functions. */
  std::vector<Matrix<float_type>> _wigner_d;

  /*! @brief Derivatives of the Wigner D functions. */
  std::vector<Matrix<float_type>> _dwigner_d;

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
   * @param precompute Precompute the Wigner D functions (this is faster, but
   * uses more memory)?
   */
  inline WignerDResources(const uint_fast32_t nmax, const uint_fast32_t ngauss,
                          const GaussBasedResources &quadrature_points,
                          const bool precompute = true)
      : _nmax(nmax), _ngauss(ngauss), _quadrature_points(quadrature_points) {

    if (precompute) {
      _wigner_d.reserve(nmax + 1);
      _dwigner_d.reserve(nmax + 1);
      for (uint_fast32_t i = 0; i < nmax + 1; ++i) {
        _wigner_d.push_back(Matrix<float_type>(2 * ngauss, nmax));
        _dwigner_d.push_back(Matrix<float_type>(2 * ngauss, nmax));
      }
    }
  }

  virtual ~WignerDResources() {}

  /**
   * @brief Get the size in memory of a hypothetical WignerDResources object
   * with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param precompute Precompute the Wigner D functions (this is faster, but
   * uses more memory)?
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ngauss,
                                       const bool precompute = true) {
    size_t size = sizeof(WignerDResources);
    if (precompute) {
      // wigner_d
      size += (nmax + 1) * 2 * ngauss * nmax * sizeof(float_type);
      // dwigner_d
      size += (nmax + 1) * 2 * ngauss * nmax * sizeof(float_type);
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
    quicksched.link_task_and_resource(*this, _quadrature_points, false);
  }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    // we only precompute the coefficients if that was allowed
    if (_wigner_d.size() > 0) {
      for (uint_fast32_t m = 0; m < _nmax + 1; ++m) {
        for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {
          const uint_fast32_t i1 = _ngauss + ig;
          const uint_fast32_t i2 = _ngauss - ig + 1;
          SpecialFunctions::wigner_dn_0m(
              _quadrature_points.get_costheta(i1 - 1), _nmax, m,
              &_wigner_d[m].get_row(i1 - 1)[0],
              &_dwigner_d[m].get_row(i1 - 1)[0]);
          int_fast8_t sign = 1;
          for (uint_fast32_t n = 0; n < _nmax; ++n) {
            sign = -sign;
            _wigner_d[m](i2 - 1, n) = sign * _wigner_d[m](i1 - 1, n);
            _dwigner_d[m](i2 - 1, n) = -sign * _dwigner_d[m](i1 - 1, n);

            ctm_assert_not_nan(_wigner_d[m](i1 - 1, n));
            ctm_assert_not_nan(_wigner_d[m](i2 - 1, n));
            ctm_assert_not_nan(_dwigner_d[m](i1 - 1, n));
            ctm_assert_not_nan(_dwigner_d[m](i2 - 1, n));
          }
        }
      }
    }
    make_available();
  }

  /**
   * @brief Get the Wigner D function for the given Gauss-Legendre quadrature
   * point.
   *
   * @param m @f$m@f$ value.
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function value.
   */
  inline float_type get_wigner_d(const uint_fast32_t m, const uint_fast32_t ig,
                                 const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();

    if (_wigner_d.size() > 0) {
      return _wigner_d[m](ig, n - 1);
    } else {
      ctm_error("Ad hoc per-element computation not available!");
    }
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given Gauss-
   * Legendre quadrature point.
   *
   * @param m @f$m@f$ value.
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param n Order, @f$n@f$.
   * @return Corresponding Wigner D function derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast32_t m, const uint_fast32_t ig,
                                  const uint_fast32_t n) const {
    // check that the resource was actually computed
    check_use();

    if (_dwigner_d.size() > 0) {
      return _dwigner_d[m](ig, n - 1);
    } else {
      ctm_error("Ad hoc per-element computation not available!");
    }
  }

  /**
   * @brief Get the Wigner D function and derivative for the given value of
   * @f$m@f$ and for the given Gauss-Legendre quadrature point.
   *
   * If the functions were precomputed, this function simply returns those
   * values. If not, the values are computed now.
   *
   * @param m @f$m@f$ value.
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @param wigner_d Pointer to an array-like structure to store the Wigner D
   * function values in (of size _nmax or more).
   * @param dwigner_d Pointer to an array-like structure to store the
   * derivatives of the Wigner D function values in (of size _nmax or more).
   */
  inline void get_wigner_functions(const uint_fast32_t m,
                                   const uint_fast32_t ig,
                                   std::vector<float_type> &wigner_d,
                                   std::vector<float_type> &dwigner_d) const {

    if (_wigner_d.size() > 0) {
      for (uint_fast32_t i = 0; i < _nmax; ++i) {
        wigner_d[i] = _wigner_d[m](ig, i);
        dwigner_d[i] = _dwigner_d[m](ig, i);
      }
    } else {
      // compute now
      SpecialFunctions::wigner_dn_0m(_quadrature_points.get_costheta(ig), _nmax,
                                     m, &wigner_d[0], &dwigner_d[0]);
    }
  }
};

#endif // WIGNERDRESOURCES_HPP
