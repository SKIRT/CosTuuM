/**
 * @file ShapeDistribution.hpp
 *
 * @brief Interface for shape distributions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SHAPEDISTRIBUTION_HPP
#define SHAPEDISTRIBUTION_HPP

#include "Configuration.hpp"

/**
 * @brief Interface for shape distributions.
 */
class ShapeDistribution {
public:
  /**
   * @brief Virtual shape distribution function.
   *
   * The default version assumes a Gaussian-like shape distribution centred
   * on an axis ratio @f$d=1@f$ (@f$d=\frac{a}{b}@f$, with @f$b@f$ the symmetry
   * axis) and limited to the range @f$[0.5, 1.5]@f$.
   *
   * @param axis_ratio Input axis ratio, @f$d = \frac{a}{b}@f$.
   * @return Value of the shape distribution function for this axis ratio.
   */
  virtual float_type operator()(const float_type axis_ratio) const {

    if (axis_ratio >= 0.5 && axis_ratio <= 1.5) {
      const float_type x = axis_ratio - 1.;
      return exp(-0.5 * x * x / 0.027777778);
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
  virtual float_type get_minimum_axis_ratio() const { return 0.5; }

  /**
   * @brief Get the maximum axis ratio @f$d@f$ for this distribution.
   *
   * The default version returns @f$1.5@f$.
   *
   * @return Maximum value for the distribution.
   */
  virtual float_type get_maximum_axis_ratio() const { return 1.5; }
};

#endif // SHAPEDISTRIBUTION_HPP
