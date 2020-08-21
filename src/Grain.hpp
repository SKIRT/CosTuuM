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
 * @file Grain.hpp
 *
 * @brief General interface for dust grains.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GRAIN_HPP
#define GRAIN_HPP

#include "Direction.hpp"
#include "IntersectionEvent.hpp"
#include "Line.hpp"
#include "Point.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief General interface for dust grains.
 */
class Grain {
public:
  /**
   * @brief Generate a random line that will intersect with the grain and that
   * travels in the given direction.
   *
   * @param direction Direction.
   * @param random_generator RandomGenerator to use.
   * @return Random line.
   */
  virtual Line
  generate_random_line(const Direction direction,
                       RandomGenerator &random_generator) const = 0;

  /**
   * @brief Get the intersection of this grain with the given line.
   *
   * @param line Line.
   * @param front Do we want the forward intersection point (default: true)?
   * @return IntersectionEvent for the intersection of this grain with the line.
   */
  virtual IntersectionEvent get_intersection(const Line line,
                                             const bool front = true) const = 0;
};

#endif // GRAIN_HPP
