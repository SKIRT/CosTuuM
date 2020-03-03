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
#include "Error.hpp"
#include "SpecialFunctions.hpp"

#include <cmath>
#include <vector>

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

protected:
  /*! @brief Points where the shape distribution is evaluated. */
  std::vector<float_type> _shapes;

  /*! @brief Weights for the evaluation points. */
  std::vector<float_type> _weights;

public:
  virtual ~ShapeDistribution() {}

  /**
   * @brief Evaluate the shape distribution at a given number of Gauss-Legendre
   * quadrature points.
   *
   * @param npoints Number of points to use.
   */
  inline void evaluate(const uint_fast32_t npoints) {

    _shapes.resize(npoints);
    _weights.resize(npoints);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        npoints, get_minimum_axis_ratio(), get_maximum_axis_ratio(), _shapes,
        _weights);

    for (uint_fast32_t i = 0; i < npoints; ++i) {
      _weights[i] *= operator()(_shapes[i]);
    }
  }

  /**
   * @brief Get the number of evaluation points for the shape distribution.
   *
   * @return Number of shapes.
   */
  inline uint_fast32_t get_number_of_points() const { return _shapes.size(); }

  /**
   * @brief Get the shape corresponding to the evaluation point with the given
   * index.
   *
   * @param index Index, needs to be smaller than get_number_of_points().
   * @return Corresponding shape value.
   */
  inline float_type get_shape(const uint_fast32_t index) const {
    ctm_assert(index < _shapes.size());
    return _shapes[index];
  }

  /**
   * @brief Get the weight corresponding to the evaluation point with the given
   * index.
   *
   * @param index Index, needs to be smaller than get_number_of_points().
   * @return Corresponding weight value.
   */
  inline float_type get_weight(const uint_fast32_t index) const {
    ctm_assert(index < _weights.size());
    return _weights[index];
  }
};

#endif // SHAPEDISTRIBUTION_HPP
