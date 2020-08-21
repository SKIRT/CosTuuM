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

#include <fstream>

/**
 * @brief Unit test for the SphericalGrain class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  SphericalGrain grain;

  /// projection test
  {

    RandomGenerator random_generator(42);

    /// "normal direction"
    {
      // projection direction
      const Direction direction(M_PI / 3., M_PI / 3.);
      // normal to the projection plane (just the reverse of the projection
      // direction)
      const Direction normal = direction.reverse();

      std::ofstream ofile("test_spherical_grain.txt");
      ofile << "# x\ty\tz\n";
      for (uint_fast32_t i = 0; i < 100000u; ++i) {
        const Line line =
            grain.generate_random_line(direction, random_generator);
        const Point base_point = line.get_base_point();
        ofile << base_point.x() << "\t" << base_point.y() << "\t"
              << base_point.z() << "\n";
        // check that the point is actually on the projection plane
        const double nr = normal.nx() * base_point.x() +
                          normal.ny() * base_point.y() +
                          normal.nz() * base_point.z();
        assert_condition(std::abs(nr) < 1.e-15);
        const Direction direction = line.get_direction();
        assert_condition(direction.get_zenith_angle() == M_PI / 3.);
        assert_condition(direction.get_azimuth_angle() == M_PI / 3.);
      }
      ofile.close();
    }

    /// direction with theta > pi/2
    {
      const Direction direction(2. * M_PI / 3., M_PI / 3.);
      // normal to the projection plane (just the reverse of the projection
      // direction)
      const Direction normal = direction.reverse();

      for (uint_fast32_t i = 0; i < 100000u; ++i) {
        const Line line =
            grain.generate_random_line(direction, random_generator);
        const Point base_point = line.get_base_point();
        // check that the point is actually on the projection plane
        const double nr = normal.nx() * base_point.x() +
                          normal.ny() * base_point.y() +
                          normal.nz() * base_point.z();
        assert_condition(std::abs(nr) < 1.e-15);
        const Direction direction = line.get_direction();
        assert_condition(direction.get_zenith_angle() == 2. * M_PI / 3.);
        assert_condition(direction.get_azimuth_angle() == M_PI / 3.);
      }
    }
  }

  /// simple intersection test
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

  /// more elaborate intersection test
  {
    RandomGenerator random_generator(32);
    for (uint_fast32_t i = 0; i < 100000u; ++i) {
      const Direction direction(
          M_PI * random_generator.get_uniform_random_double(),
          -M_PI + 2. * M_PI * random_generator.get_uniform_random_double());
      ctm_warning("direction: %g %g", direction.get_zenith_angle(),
                  direction.get_azimuth_angle());
      const Line line = grain.generate_random_line(direction, random_generator);
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
