/**
 * @file OrientationDistribution.hpp
 *
 * @brief Class that contains the coefficients of the Legendre function
 * expansion of an orientation distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ORIENTATIONDISTRIBUTION_HPP
#define ORIENTATIONDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "SpecialFunctions.hpp"

/**
 * @brief Class that contains the coefficients of the Legendre function
 * expansion of an orientation distribution.
 */
class OrientationDistribution {
private:
  /*! @brief Highest order for which we store a coefficient. */
  const uint_fast32_t _maximum_order;

  /*! @brief Coefficients. */
  std::vector<float_type> _coefficients;

  /**
   * @brief Integrand for coefficients.
   *
   * @param beta Angle, @f$\beta{}@f$.
   * @param args Additional arguments: the order of the Wigner D function in the
   * integrand.
   * @return Value of the integrand.
   */
  static float_type integrand(const float_type beta, void *args) {

    const uint_fast32_t order = *reinterpret_cast<uint_fast32_t *>(args);

    const float_type cosbeta = cos(beta);
    const float_type sinbeta = sin(beta);
    float_type wigner_d_n;
    if (order > 0) {
      std::vector<float_type> wigner_d(order), dwigner_d(order);
      SpecialFunctions::wigner_dn_0m(cosbeta, order, 0, &wigner_d[0],
                                     &dwigner_d[0]);
      wigner_d_n = wigner_d[order - 1];
    } else {
      wigner_d_n = 1.;
    }
    const float_type px = 0.5 + 0.1 * (cosbeta * cosbeta - 1.);
    return sinbeta * px * wigner_d_n;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   */
  OrientationDistribution(const uint_fast32_t nmax)
      : _maximum_order(nmax + 1), _coefficients(_maximum_order, 0.) {

    for (uint_fast32_t i = 0; i < _coefficients.size(); ++i) {
      uint_fast32_t order = i;
      ctm_warning("n = %" PRIuFAST32, order);
      _coefficients[i] =
          SpecialFunctions::gauss_legendre_quadrature<float_type>(
              integrand, 0., M_PI, &order, 2 * (order + 1), 100 * (order + 1),
              1.e-10, 1.e-5);
    }
  }

  /**
   * @brief Get coefficient @f$p_n@f$.
   *
   * @param n Order @f$n@f$.
   * @return Corresponding coefficient.
   */
  inline float_type get_coefficient(const uint_fast32_t n) const {
    return _coefficients[n];
  }
};

#endif // ORIENTATIONDISTRIBUTION_HPP
