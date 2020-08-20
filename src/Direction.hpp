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
};

#endif // DIRECTION_HPP
