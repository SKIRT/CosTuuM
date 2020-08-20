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
 * @file Point.hpp
 *
 * @brief Representation of a 3D point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef POINT_HPP
#define POINT_HPP

#include "Direction.hpp"

#include <cmath>

/**
 * @brief Representation of a 3D point.
 */
class Point {
private:
  /*! @brief Coordinates of the point. */
  double _x[3];

public:
  /**
   * @brief Constructor.
   *
   * @param x X coordinate.
   * @param y Y coordinate.
   * @param z Z coordinate.
   */
  inline Point(const double x, const double y, const double z) : _x{x, y, z} {}

  /**
   * @brief Get the x-coordinate of the point.
   *
   * @return X-coordinate.
   */
  inline double x() const { return _x[0]; }

  /**
   * @brief Get the y-coordinate of the point.
   *
   * @return Y-coordinate.
   */
  inline double y() const { return _x[1]; }

  /**
   * @brief Get the z-coordinate of the point.
   *
   * @return Z-coordinate.
   */
  inline double z() const { return _x[2]; }

  /**
   * @brief Convert the coordinates of the point to spherical coordinates.
   *
   * @param r Output radius.
   * @param theta Output zenith angle.
   * @param phi Output azimuth angle.
   */
  inline void spherical_coordinates(double &r, double &theta,
                                    double &phi) const {
    r = std::sqrt(_x[0] * _x[0] + _x[1] * _x[1] + _x[2] * _x[2]);
    if (r != 0.) {
      theta = std::acos(_x[2] / r);
      phi = std::atan2(_x[1], _x[0]);
    } else {
      theta = 0.;
      phi = 0.;
    }
  }

  /**
   * @brief Return a new Point that is obtained by translating this point over
   * the given distance along the given direction.
   *
   * @param t Distance.
   * @param direction Direction.
   * @return New point.
   */
  inline Point translate(const double t, const Direction direction) const {
    return Point(_x[0] + t * direction.nx(), _x[1] + t * direction.ny(),
                 _x[2] + t * direction.nz());
  }
};

#endif // POINT_HPP
