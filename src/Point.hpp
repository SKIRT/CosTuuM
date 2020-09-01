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

  /**
   * @brief Scale the point coordinates with the given scale factor.
   *
   * @param scale_factor Scale factor.
   * @return Scaled point.
   */
  inline Point scale(const double scale_factor) const {
    return Point(scale_factor * _x[0], scale_factor * _x[1],
                 scale_factor * _x[2]);
  }

  /**
   * @brief Get the cross product of the point vector with the given direction.
   *
   * @param direction Direction.
   * @return Cross product.
   */
  inline Point cross(const Direction direction) const {

    const double x = _x[1] * direction.nz() - _x[2] * direction.ny();
    const double y = _x[2] * direction.nx() - _x[0] * direction.nz();
    const double z = _x[0] * direction.ny() - _x[1] * direction.nx();

    return Point(x, y, z);
  }

  /**
   * @brief Get the dot product of the point vector and the given direction.
   *
   * @param direction Direction.
   * @return Dot product.
   */
  inline double dot(const Direction direction) const {
    return _x[0] * direction.nx() + _x[1] * direction.ny() +
           _x[2] * direction.nz();
  }

  /**
   * @brief Rotate the point vector over the given angle around the given axis,
   * according to the right hand rule.
   *
   * It is assumed that the axis is perpendicular to the point vector.
   *
   * @param axis Rotation axis.
   * @param angle Rotation angle.
   * @return Rotated direction.
   */
  inline Point rotate(const Direction axis, const double angle) const {

    const Point p_cross_a = cross(axis);
    const double pdota = dot(axis);
    const double cosa = std::cos(angle);
    const double sina = std::sin(angle);
    const double factor = pdota * (1. - cosa);

    const double x = _x[0] * cosa - p_cross_a._x[0] * sina + axis.nx() * factor;
    const double y = _x[1] * cosa - p_cross_a._x[1] * sina + axis.ny() * factor;
    const double z = _x[2] * cosa - p_cross_a._x[2] * sina + axis.nz() * factor;

    return Point(x, y, z);
  }
};

#endif // POINT_HPP
