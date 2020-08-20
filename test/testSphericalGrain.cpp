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
 * @file testSphericalGrain.cpp
 *
 * @brief Unit test for the SphericalGrain class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "SphericalGrain.hpp"

/**
 * @brief Unit test for the SphericalGrain class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  SphericalGrain grain;

  /// known, simple test case
  {
    const Point point(0.5, 0.5, 0.);
    const Direction direction(M_PI, 0.);
    const Line line(point, direction);

    // backwards
    {
      const IntersectionEvent event = grain.get_intersection(line, false);
      const Point intersection_point = event.get_intersection_point();
      ctm_warning("intersection: %g %g %g", intersection_point.x(),
                  intersection_point.y(), intersection_point.z());
      assert_values_equal_rel(intersection_point.x(), 0.5, 1.e-15);
      assert_values_equal_rel(intersection_point.y(), 0.5, 1.e-15);
      assert_values_equal_rel(intersection_point.z(), std::sqrt(0.5), 1.e-15);
    }
    // front
    {
      const IntersectionEvent event = grain.get_intersection(line, true);
      const Point intersection_point = event.get_intersection_point();
      ctm_warning("intersection: %g %g %g", intersection_point.x(),
                  intersection_point.y(), intersection_point.z());
      assert_values_equal_rel(intersection_point.x(), 0.5, 1.e-15);
      assert_values_equal_rel(intersection_point.y(), 0.5, 1.e-15);
      assert_values_equal_rel(intersection_point.z(), -std::sqrt(0.5), 1.e-15);
    }
  }

  /// more elaborate test
  {
    RandomGenerator random_generator(32);
    for (uint_fast32_t i = 0; i < 100000u; ++i) {
      const Direction direction(
          M_PI * random_generator.get_uniform_random_double(),
          -M_PI + 2. * M_PI * random_generator.get_uniform_random_double());
      ctm_warning("direction: %g %g", direction.get_zenith_angle(),
                  direction.get_azimuth_angle());
      const GrainProjection *projection = grain.get_projection(direction);
      const Line line = projection->generate_random_line(random_generator);
      const IntersectionEvent event = grain.get_intersection(line, true);
      const Point intersection_point = event.get_intersection_point();
      ctm_warning("intersection: %g %g %g", intersection_point.x(),
                  intersection_point.y(), intersection_point.z());
      double r, theta, phi;
      intersection_point.spherical_coordinates(r, theta, phi);
      ctm_warning("intersection: %g %g %g", r, theta, phi);
      assert_values_equal_rel(r, 1., 1.e-15);
      assert_condition(std::abs(theta - direction.get_zenith_angle()) <=
                       0.5 * M_PI);
    }
  }

  return 0;
}
