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
class WignerDResources : public Resource, public Task {
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
  GaussBasedResources &_quadrature_points;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param quadrature_points Gauss-Legendre quadrature points (to be computed
   * before this task is executed!).
   */
  inline WignerDResources(const uint_fast32_t nmax, const uint_fast32_t ngauss,
                          GaussBasedResources &quadrature_points)
      : _nmax(nmax), _ngauss(ngauss), _wigner_d(2 * ngauss, nmax),
        _dwigner_d(2 * ngauss, nmax), _quadrature_points(quadrature_points) {}

  virtual ~WignerDResources() {}

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
      }
    }
  }
};

#endif // WIGNERDRESOURCES_HPP
