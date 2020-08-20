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
 * @file IntersectionEvent.hpp
 *
 * @brief Object storing information about the intersection of a Grain with a
 * Line.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef INTERSECTIONEVENT_HPP
#define INTERSECTIONEVENT_HPP

#include "Direction.hpp"
#include "Point.hpp"

/**
 * @brief Object storing information about the intersection of a Grain with a
 * Line.
 */
class IntersectionEvent {
private:
  /*! @brief Intersection point. */
  const Point _point;

  /*! @brief Normal at the intersection point. */
  const Direction _normal;

public:
  /**
   * @brief Constructor.
   *
   * @param point Intersection point.
   * @param normal Normal at the intersection point.
   */
  inline IntersectionEvent(const Point point, const Direction normal)
      : _point(point), _normal(normal) {}

  /**
   * @brief Get the intersection point.
   *
   * @return Intersection point.
   */
  inline Point get_intersection_point() const { return _point; }

  /**
   * @brief Get the normal at the intersection point.
   *
   * @return Normal.
   */
  inline Direction get_normal() const { return _normal; }
};

#endif // INTERSECTIONEVENT_HPP
