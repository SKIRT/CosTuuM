/**
 * @file DraineHensleyShapeDistribution.hpp
 *
 * @brief Draine & Hensley (2017) CDE2 shape distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DRAINEHENSLEYSHAPEDISTRIBUTION_HPP
#define DRAINEHENSLEYSHAPEDISTRIBUTION_HPP

#include "Error.hpp"
#include "ShapeDistribution.hpp"

/**
 * @brief Draine & Hensley (2017) CDE2 shape distribution.
 */
class DraineHensleyShapeDistribution : public ShapeDistribution {
private:
  /**
   * @brief Get the shape factor @f$L@f$ corresponding to the given axis ratio
   * @f$d@f$.
   *
   * We use the expressions from Min et al. (2013)
   * (https://ui.adsabs.harvard.edu/abs/2003A%26A...404...35M/abstract), eq.
   * 22:
   * @f[
   *   L = \frac{1-e^2}{e^2} \left[
         \frac{1}{2e} \ln\left(\frac{1+e}{1-e}\right) -1 \right],
   * @f]
   * @f[
   *   e^2 = 1 - d^2
   * @f]
   * for prolate spheroids (@f$d<1@f$),
   * @f[
   *   L = \frac{1}{e^2} \left[
   *     1 - \frac{\sqrt{1-e^2}}{e}\arcsin(e) \right],
   * @f]
   * @f[
   *   e^2 = 1 - \frac{1}{d^2}
   * @f]
   * for oblate spheroids (@f$d>1@f$).
   *
   * @param axis_ratio Axis ratio, @f$d@f$.
   * @return Corresponding shape factor, @f$L(d)@f$.
   */
  inline static float_type get_shape_factor(const float_type axis_ratio) {
    if (axis_ratio < 1.) {
      // prolate spheroid
      const float_type e2 = 1. - axis_ratio * axis_ratio;
      const float_type e = sqrt(e2);
      return (1. - e2) * (0.5 * log((1. + e) / (1. - e)) / e - 1.) / e2;
    } else if (axis_ratio > 1.) {
      // oblate spheroid
      const float_type e2 = 1. - 1. / (axis_ratio * axis_ratio);
      const float_type e = sqrt(e2);
      return (1. - sqrt(1. - e2) * asin(e) / e) / e2;
    } else {
      ctm_assert(axis_ratio == 1.);
      // sphere
      return 1. / 3.;
    }
  }

  /**
   * @brief Get the Jacobian transformation @f$\frac{{\rm{}d}L}{{\rm{}d}d}@f$
   * for the given axis ratio @f$d@f$.
   *
   * We use the expressions from Min et al. (2013)
   * (https://ui.adsabs.harvard.edu/abs/2003A%26A...404...35M/abstract), eq.
   * 25:
   * @f[
   *   \frac{{\rm{}d}L}{{\rm{}d}d} = \frac{\sqrt{1-e^2}}{2e^5} \left[
   *      (3-e^2) \ln\left(\frac{1+e}{1-e}\right) -6e \right],
   * @f]
   * @f[
   *   e^2 = 1 - d^2
   * @f]
   * for prolate spheroids (@f$d<1@f$),
   * @f[
   *   \frac{{\rm{}d}L}{{\rm{}d}d} = \frac{1-e^2}{e^5} \left[
   *     (3-2e^2) \arcsin(e) - 3e\sqrt{1-e^2} \right],
   * @f]
   * @f[
   *   e^2 = 1 - \frac{1}{d^2}
   * @f]
   * for oblate spheroids (@f$d>1@f$).
   *
   * These expressions are divergent for @f$e\rightarrow{}0@f$, although it can
   * be shown from the Taylor expansions around @f$e=0@f$:
   * @f[
   *   \frac{{\rm{}d}L}{{\rm{}d}d} = \frac{4}{15} + \frac{2}{21} e^2
   *     + \frac{3}{70} e^4 + \dots{}
   * @f]
   * (for prolate spheroids)
   * @f[
   *   \frac{{\rm{}d}L}{{\rm{}d}d} = \frac{4}{15} - \frac{2}{21} e^2
   *     - \frac{11}{210} e^4 + \dots{}
   * @f]
   * that both have the same limit. We will use these series expansions to deal
   * with value of @f$e@f$ close to @f$0@f$.
   *
   * @param axis_ratio Axis ratio, @f$d@f$.
   * @return Jacobian, @f$\frac{{\rm{}d}L}{{\rm{}d}d}(d)@f$.
   */
  inline static float_type get_jacobian(const float_type axis_ratio) {
    if (axis_ratio < 1.) {
      // prolate spheroid
      const float_type e2 = 1. - axis_ratio * axis_ratio;
      if (e2 < 1.e-5) {
        // Taylor expansion
        return 4. / 15. + 2. * e2 / 21. + 3. * e2 * e2 / 70.;
      } else {
        // full expression
        const float_type e = sqrt(e2);
        return 0.5 * sqrt(1. - e2) *
               ((3. - e2) * log((1. + e) / (1. - e)) - 6. * e) / (e2 * e2 * e);
      }
    } else if (axis_ratio > 1.) {
      // oblate spheroid
      const float_type e2 = 1. - 1. / (axis_ratio * axis_ratio);
      if (e2 < 1.e-5) {
        // Taylor expansion
        return 4. / 15. - 2. * e2 / 21. - 11. * e2 * e2 / 210.;
      } else {
        // full expression
        const float_type e = sqrt(e2);
        return (1. - e2) * ((3. - 2. * e2) * asin(e) - 3. * e * sqrt(1. - e2)) /
               (e2 * e2 * e);
      }
    } else {
      ctm_assert(axis_ratio == 1.);
      // sphere
      return 4. / 15.;
    }
  }

  /**
   * @brief Get the value of the CDE2 distribution for the given shape factor.
   *
   * The general expression for the CDE2 distribution is given in Draine &
   * Hensley (2017)
   * (https://ui.adsabs.harvard.edu/abs/2017arXiv171008968D/abstract), eq. 22:
   * @f[
   *   G(L_1, L_2) = 120 L_1 L_2 L_3 = 120 L_1 L_2 (1 - L_1 - L_2).
   * @f]
   * In this expression, the shape factors @f$L_1 \geq{} L_2 \geq{} L_3@f$
   * (@f$L_1+L_2+L_3=1@f$) are generally not the same. For spheroids, two of
   * the shape factors will be the same, so we can reparametrise the
   * distribution as
   * @f[
   *   G(L) = 12 L (1-L)^2,
   * @f]
   * where the normalisation was chosen so that
   * @f$\int_0^1 G(L) {\rm{}d}L = 1@f$.
   *
   * In this reparametrisation, @f$L@f$ takes respectively the role of the
   * long axis shape factor @f$L_3@f$ for prolate spheroids
   * (@f$L<\frac{1}{3}@f$)and the role of the short axis shape factor @f$L_1@f$
   * for oblate spheroids (@f$L>\frac{1}{3}@f$).
   *
   * @param L Shape factor, @f$L@f$.
   * @return Corresponding CDE2 distribution function, @f$G(L)@f$.
   */
  inline static float_type cde2(const float_type L) {
    const float_type Lm1 = L - 1.;
    return 12. * L * Lm1 * Lm1;
  }

public:
  /**
   * @brief Virtual shape distribution function.
   *
   * @param axis_ratio Input axis ratio, @f$d = \frac{a}{b}@f$.
   * @return Value of the shape distribution function for this axis ratio.
   */
  virtual float_type operator()(const float_type axis_ratio) const {

    if (axis_ratio >= 0.12 && axis_ratio <= 3.68) {
      return cde2(get_shape_factor(axis_ratio)) * get_jacobian(axis_ratio);
    } else {
      return 0.;
    }
  }

  /**
   * @brief Get the minimum axis ratio @f$d@f$ for this distribution.
   *
   * The default version returns @f$0.5@f$.
   *
   * @return Minimum value for the distribution.
   */
  virtual float_type get_minimum_axis_ratio() const { return 0.12; }

  /**
   * @brief Get the maximum axis ratio @f$d@f$ for this distribution.
   *
   * The default version returns @f$1.5@f$.
   *
   * @return Maximum value for the distribution.
   */
  virtual float_type get_maximum_axis_ratio() const { return 3.68; }
};

#endif // DRAINEHENSLEYSHAPEDISTRIBUTION_HPP
