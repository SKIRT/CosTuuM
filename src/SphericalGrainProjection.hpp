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
 * @file SphericalGrainProjection.hpp
 *
 * @brief GrainProjection implementation for spherical grains.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPHERICALGRAINPROJECTION_HPP
#define SPHERICALGRAINPROJECTION_HPP

#include "GrainProjection.hpp"

#include <cmath>

/**
 * @brief GrainProjection implementation for spherical grains.
 */
class SphericalGrainProjection : public GrainProjection {
public:
  /**
   * @brief Constructor.
   *
   * @param direction Direction.
   */
  inline SphericalGrainProjection(const Direction direction)
      : GrainProjection(direction) {}

  /**
   * @brief Generate a random line that will intersect with the grain.
   *
   * @param random_generator RandomGenerator to use.
   * @return Random line.
   */
  virtual Line generate_random_line(RandomGenerator &random_generator) const {

    // first generate a random radius and polar angle in the plane perpendicular
    // to the projection direction
    const double r = std::sqrt(random_generator.get_uniform_random_double());
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();

    // now convert these to spherical coordinates
    // the radius (distance between sphere/circle origin and random point) stays
    // the same
    // the zenith angle is fixed by the requirement that the plane is
    // perpendicular to the projection direction
    // the azimuth angle can be chosen arbitrary in this plane, since we have
    // azimuthal symmetry

    // compute the fixed zenith angle
    double theta = 0.5 * M_PI - _direction.get_zenith_angle();
    // make sure theta \in [0, pi[
    // note that each of these theta flips would in principle require a
    // corresponding rotation of the azimuth angle over an angle pi
    // however, due to azimuthal symmetry, we do not actually need to perform
    // these rotations
    if (theta < 0.) {
      theta = -theta;
    }
    if (theta > M_PI) {
      theta = 2. * M_PI - theta;
    }

    // now compute the cartesian coordinates
    const double sintheta = std::sin(theta);
    const double x = r * sintheta * std::cos(phi);
    const double y = r * sintheta * std::sin(phi);
    const double z = r * std::cos(theta);

    // construct the line
    const Point point(x, y, z);
    return Line(point, _direction);
  }
};

#endif // SPHERICALGRAINPROJECTION_HPP
