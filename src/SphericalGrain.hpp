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
 * @file SphericalGrain.hpp
 *
 * @brief Grain implementation for spherical grains.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPHERICALGRAIN_HPP
#define SPHERICALGRAIN_HPP

#include "Error.hpp"
#include "Grain.hpp"
#include "SphericalGrainProjection.hpp"

/**
 * @brief Grain implementation for spherical grains.
 */
class SphericalGrain : public Grain {
public:
  /**
   * @brief Get the projection of this grain perpendicular to the given
   * direction.
   *
   * @param direction Direction.
   * @return Projection perpendicular to this direction.
   */
  virtual GrainProjection *get_projection(const Direction direction) const {
    return new SphericalGrainProjection(direction);
  }

  /**
   * @brief Get the intersection of this grain with the given line.
   *
   * @param line Line.
   * @param front Do we want the forward intersection point (default: true)?
   * @return IntersectionEvent for the intersection of this grain with the line.
   */
  virtual IntersectionEvent get_intersection(const Line line,
                                             const bool front) const {

    const Point p0 = line.get_base_point();
    const Direction d = line.get_direction();
    const double xnx = p0.x() * d.nx();
    const double yny = p0.y() * d.ny();
    const double znz = p0.z() * d.nz();
    const double discriminant = xnx * xnx + yny * yny + znz * znz +
                                2. * (xnx * yny + xnx * znz + yny * znz) -
                                p0.x() * p0.x() - p0.y() * p0.y() -
                                p0.z() * p0.z() + 1.;
    double t;
    if (discriminant < 0.) {
      ctm_error("No valid intersection point!");
      t = 0.;
    } else if (discriminant == 0.) {
      // only 1 intersection point
      t = -xnx - yny - znz;
    } else {
      if (front) {
        t = std::sqrt(discriminant) - xnx - yny - znz;
      } else {
        t = -std::sqrt(discriminant) - xnx - yny - znz;
      }
    }

    const Point intersection_point = line.evaluate(t);

    // the normal in the intersection point happens to be the direction of the
    // intersection point (since this is already a unit sphere)
    double r, theta, phi;
    intersection_point.spherical_coordinates(r, theta, phi);
    const Direction normal(theta, phi);

    return IntersectionEvent(intersection_point, normal);
  }
};

#endif // SPHERICALGRAIN_HPP
