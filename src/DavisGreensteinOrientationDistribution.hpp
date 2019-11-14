/**
 * @file DavisGreensteinOrientationDistribution.hpp
 *
 * @brief Mishchenko (1991) perfect Davis-Greenstein alignment orientation
 * distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DAVISGREENSTEINORIENTATIONDISTRIBUTION_HPP
#define DAVISGREENSTEINORIENTATIONDISTRIBUTION_HPP

#include "OrientationDistribution.hpp"

/**
 * @brief Mishchenko (1991) perfect Davis-Greenstein alignment orientation
 * distribution.
 *
 * Based on Mishchenko, 1991, ApJ, 367, 561 (https://doi.org/10.1086/169652).
 *
 * This is a trivial distribution function for which no numerical integration
 * for the coefficients is required.
 *
 * For oblate spheroidals, the distribution is given by
 * @f[
 *  p(\beta{}) = \delta(\cos(\beta{}) - 1),
 * @f]
 * with @f$\delta(x)@f$ the Dirac-delta distribution. The corresponding
 * expansion coefficients are
 * @f[
 *  p_n = 1, \forall{} n.
 * @f]
 * For prolate spheroidals the distribution function is
 * @f[
 *  p(\beta{}) = \delta(\cos(\beta{})),
 * @f]
 * with odd coefficients
 * @f[
 *  p_{2n+1} = 0, \forall{} n
 * @f]
 * and even coefficients
 * @f[
 *  p_{2n} = (-1)^n \frac{1\dot{}3\dotso{}(2n-1)}{2\dot{}4\dotso{}2n},
 * \forall{} n.
 * @f]
 * Note that there is some ambiguity for the even component and @f$n=0@f$.
 * Setting @f$p_0 = 1@f$ reproduces the results in Table 1 of Mishchenko (1991).
 */
class DavisGreensteinOrientationDistribution : public OrientationDistribution {
public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   * @param axis_ratio Axis ratio @f$d@f$ of the spheroidal scatterers, used to
   * determine whether the prolate (@f$d < 1@f$) or oblate (@f$d > 1@f$)
   * parameter formulae should be used.
   */
  inline DavisGreensteinOrientationDistribution(const uint_fast32_t nmax,
                                                const float_type axis_ratio)
      : OrientationDistribution(nmax) {

    if (axis_ratio < 1.) {
      // prolate spheroid, particle aligned perpendicular to the field
      // first do the odd components...
      for (uint_fast32_t i = 1; i < _coefficients.size(); i += 2) {
        _coefficients[i] = 0.;
      }
      // ...now the even components
      float_type nominator = 1.;
      float_type denominator = 2.;
      _coefficients[0] = 1.;
      for (uint_fast32_t i = 1; 2 * i < _coefficients.size(); ++i) {
        const int_fast8_t sign = (i % 2 == 0) ? 1 : -1;
        _coefficients[2 * i] = sign * nominator / denominator;
        nominator *= 2. * i + 1.;
        denominator *= 2. * (i + 1.);
      }
    } else {
      // oblate spheroid, particle aligned with the field
      for (uint_fast32_t i = 0; i < _coefficients.size(); ++i) {
        _coefficients[i] = 1.;
      }
    }
  }

  virtual ~DavisGreensteinOrientationDistribution() {}
};

#endif // DAVISGREENSTEINORIENTATIONDISTRIBUTION_HPP
