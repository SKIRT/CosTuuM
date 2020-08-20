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

/**
 * @brief 3D direction vector.
 */
class Direction {
private:
  /*! @brief Zenith angle. */
  double _theta;

  /*! @brief Azimuth angle. */
  double _phi;

public:
  /**
   * @brief Constructor.
   *
   * @param theta Zenith angle.
   * @param phi Azimuth angle.
   */
  inline Direction(const double theta, const double phi)
      : _theta(theta), _phi(phi) {}

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
};

#endif // DIRECTION_HPP
