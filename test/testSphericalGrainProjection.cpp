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
 * @file testSphericalGrainProjection.cpp
 *
 * @brief Unit test for the SphericalGrainProjection class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "SphericalGrainProjection.hpp"

#include <fstream>

/**
 * @brief Unit test for the SphericalGrainProjection class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  RandomGenerator random_generator(42);

  /// "normal direction"
  {
    const Direction direction(M_PI / 3., M_PI / 3.);
    const SphericalGrainProjection projection(direction);

    std::ofstream ofile("test_spherical_grain_projection.txt");
    ofile << "# x\ty\tz\n";
    for (uint_fast32_t i = 0; i < 100000u; ++i) {
      const Line line = projection.generate_random_line(random_generator);
      const Point base_point = line.get_base_point();
      ofile << base_point.x() << "\t" << base_point.y() << "\t"
            << base_point.z() << "\n";
      double r, theta, phi;
      base_point.spherical_coordinates(r, theta, phi);
      assert_condition(r <= 1.);
      if (r != 0.) {
        assert_values_equal_rel(theta, M_PI / 6., 1.e-15);
      }
      const Direction direction = line.get_direction();
      assert_condition(direction.get_zenith_angle() == M_PI / 3.);
      assert_condition(direction.get_azimuth_angle() == M_PI / 3.);
    }
    ofile.close();
  }

  /// direction with theta > pi/2
  {
    const Direction direction(2. * M_PI / 3., M_PI / 3.);
    const SphericalGrainProjection projection(direction);

    for (uint_fast32_t i = 0; i < 100000u; ++i) {
      const Line line = projection.generate_random_line(random_generator);
      const Point base_point = line.get_base_point();
      double r, theta, phi;
      base_point.spherical_coordinates(r, theta, phi);
      assert_condition(r <= 1.);
      if (r != 0.) {
        assert_values_equal_rel(theta, M_PI / 6., 1.e-15);
      }
      const Direction direction = line.get_direction();
      assert_condition(direction.get_zenith_angle() == 2. * M_PI / 3.);
      assert_condition(direction.get_azimuth_angle() == M_PI / 3.);
    }
  }

  return 0;
}
