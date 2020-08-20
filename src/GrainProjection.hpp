/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file GrainProjection.hpp
 *
 * @brief General interface for dust grain projections.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GRAINPROJECTION_HPP
#define GRAINPROJECTION_HPP

#include "Direction.hpp"
#include "Line.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief General interface for dust grain projections.
 */
class GrainProjection {
protected:
  /*! @brief Viewing Direction for this projection. */
  const Direction _direction;

public:
  /**
   * @brief Constructor.
   *
   * @param direction Viewing Direction for this projection.
   */
  inline GrainProjection(const Direction direction) : _direction(direction) {}

  /**
   * @brief Generate a random line that will intersect with the grain.
   *
   * @param random_generator RandomGenerator to use.
   * @return Random line.
   */
  virtual Line
  generate_random_line(RandomGenerator &random_generator) const = 0;
};

#endif // GRAINPROJECTION_HPP
