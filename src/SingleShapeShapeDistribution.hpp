/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file SingleShapeShapeDistribution.hpp
 *
 * @brief Shape distribution that returns a single shape.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SINGLESHAPESHAPEDISTRIBUTION_HPP
#define SINGLESHAPESHAPEDISTRIBUTION_HPP

#include "ShapeDistribution.hpp"

/**
 * @brief Interface for shape distributions.
 */
class SingleShapeShapeDistribution : public ShapeDistribution {
public:
  /**
   * @brief Virtual shape distribution function.
   *
   * @param axis_ratio Input axis ratio, @f$d = \frac{a}{b}@f$.
   * @return Value of the shape distribution function for this axis ratio.
   */
  virtual float_type operator()(const float_type axis_ratio) const {

    if (axis_ratio == _shapes[0]) {
      return 1.;
    } else {
      return 0.;
    }
  }

  /**
   * @brief Get the minimum axis ratio @f$d@f$ for this distribution.
   *
   * @return Minimum value for the distribution.
   */
  virtual float_type get_minimum_axis_ratio() const { return _shapes[0]; }

  /**
   * @brief Get the maximum axis ratio @f$d@f$ for this distribution.
   *
   * @return Maximum value for the distribution.
   */
  virtual float_type get_maximum_axis_ratio() const { return _shapes[0]; }

  /**
   * @brief Constructor.
   *
   * @param axis_ratio Single axis ratio value returned by this distribution.
   */
  inline SingleShapeShapeDistribution(const float_type axis_ratio) {
    _shapes.resize(1, axis_ratio);
    _weights.resize(1, 1.);
  }

  virtual ~SingleShapeShapeDistribution() {}
};

#endif // SINGLESHAPESHAPEDISTRIBUTION_HPP
