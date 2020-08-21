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
 * @file Direction.hpp
 *
 * @brief 3D direction vector.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DIRECTION_HPP
#define DIRECTION_HPP

#include <cmath>

/**
 * @brief 3D direction vector.
 */
class Direction {
private:
  /*! @brief Zenith angle. */
  double _theta;

  /*! @brief Azimuth angle. */
  double _phi;

  /*! @brief Direction vector. */
  double _n[3];

public:
  /**
   * @brief Constructor.
   *
   * @param theta Zenith angle.
   * @param phi Azimuth angle.
   */
  inline Direction(const double theta, const double phi)
      : _theta(theta), _phi(phi) {

    const double sintheta = std::sin(theta);
    _n[0] = sintheta * std::cos(phi);
    _n[1] = sintheta * std::sin(phi);
    _n[2] = std::cos(theta);
  }

  /**
   * @brief Constructor.
   *
   * @param nx X-component.
   * @param ny Y-component.
   * @param nz Z-component.
   */
  inline Direction(const double nx, const double ny, const double nz)
      : _n{nx, ny, nz} {

    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    _n[0] /= norm;
    _n[1] /= norm;
    _n[2] /= norm;

    _theta = std::acos(_n[2]);
    _phi = std::atan2(_n[1], _n[0]);
  }

  /**
   * @brief Get the zenith angle.
   *
   * @return Zenith angle.
   */
  inline double get_zenith_angle() const { return _theta; }

  /**
   * @brief Get the azimuth angle.
   *
   * @return Azimuth angle.
   */
  inline double get_azimuth_angle() const { return _phi; }

  /**
   * @brief Get the x-component of the direction vector.
   *
   * @return X-component of the direction vector.
   */
  inline double nx() const { return _n[0]; }

  /**
   * @brief Get the y-component of the direction vector.
   *
   * @return Y-component of the direction vector.
   */
  inline double ny() const { return _n[1]; }

  /**
   * @brief Get the z-component of the direction vector.
   *
   * @return Z-component of the direction vector.
   */
  inline double nz() const { return _n[2]; }

  /**
   * @brief Get the angle between this direction and the given direction.
   *
   * @param direction Other direction.
   * @return Angle between the two directions.
   */
  inline double angle(const Direction direction) const {
    return std::acos(_n[0] * direction._n[0] + _n[1] * direction._n[1] +
                     _n[2] * direction._n[2]);
  }

  /**
   * @brief Get the cross product of this direction with the given direction.
   *
   * @param direction Other direction.
   * @return Cross product.
   */
  inline Direction cross(const Direction direction) const {

    const double nx = _n[1] * direction._n[2] - _n[2] * direction._n[1];
    const double ny = _n[2] * direction._n[0] - _n[0] * direction._n[2];
    const double nz = _n[0] * direction._n[1] - _n[1] * direction._n[0];

    return Direction(nx, ny, nz);
  }

  /**
   * @brief Return the reverse direction.
   *
   * @return Reverse direction.
   */
  inline Direction reverse() const { return Direction(-_n[0], -_n[1], -_n[2]); }

  /**
   * @brief Rotate the direction over the given angle around the given axis,
   * according to the right hand rule.
   *
   * It is assumed that the axis is perpendicular to the direction.
   *
   * @param axis Rotation axis.
   * @param angle Rotation angle.
   * @return Rotated direction.
   */
  inline Direction rotate_perpendicular(const Direction axis,
                                        const double angle) const {

    const Direction d_cross_s = cross(axis);

    const double cosa = std::cos(angle);
    const double sina = std::sin(angle);

    const double nx = _n[0] * cosa - d_cross_s._n[0] * sina;
    const double ny = _n[1] * cosa - d_cross_s._n[1] * sina;
    const double nz = _n[2] * cosa - d_cross_s._n[2] * sina;

    return Direction(nx, ny, nz);
  }
};

#endif // DIRECTION_HPP
