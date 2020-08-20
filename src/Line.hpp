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
 * @file Line.hpp
 *
 * @brief Representation of a 3D line.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LINE_HPP
#define LINE_HPP

#include "Direction.hpp"
#include "Point.hpp"

/**
 * @brief Representation of a 3D line.
 */
class Line {
private:
  /*! @brief Representative point on the line. */
  Point _point;

  /*! @brief Direction of the line. */
  Direction _direction;

public:
  /**
   * @brief Constructor.
   *
   * @param point Representative point on the line.
   * @param direction Direction of the line.
   */
  inline Line(const Point point, const Direction direction)
      : _point(point), _direction(direction) {}

  /**
   * @brief Get the base point for this line.
   *
   * @return Base point (point that is stored internally).
   */
  inline Point get_base_point() const { return _point; }

  /**
   * @brief Get the direction of the line.
   *
   * @return Direction (that is stored internally).
   */
  inline Direction get_direction() const { return _direction; }
};

#endif // LINE_HPP
