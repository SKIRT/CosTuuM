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
 * @file testPoint.cpp
 *
 * @brief Unit test for the Point class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "Error.hpp"
#include "Point.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Unit test for the Point class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  /// normal point
  {
    const Point point(1., 1., 1.);

    double r, theta, phi;
    point.spherical_coordinates(r, theta, phi);
    assert_values_equal_rel(r, std::sqrt(3.), 1.e-19);
    assert_values_equal_rel(theta, std::acos(1. / std::sqrt(3.)), 1.e-19);
    assert_values_equal_rel(phi, 0.25 * M_PI, 1.e-19);

    const Direction direction(0., 0., 1.);
    const Point translated_point = point.translate(1., direction);
    assert_condition(translated_point.x() == point.x());
    assert_condition(translated_point.y() == point.y());
    assert_condition(translated_point.z() == point.z() + 1.);

    const Point rotated_point = point.rotate(direction, 0.5 * M_PI);
    ctm_warning("rotated point: %g %g %g", rotated_point.x(), rotated_point.y(),
                rotated_point.z());
    assert_values_equal_rel(rotated_point.x(), -1., 1.e-16);
    assert_values_equal_rel(rotated_point.y(), 1., 1.e-16);
    assert_condition(rotated_point.z() == point.z());
  }

  /// special point
  {
    const Point point(-1., -1., -1.);

    double r, theta, phi;
    point.spherical_coordinates(r, theta, phi);
    assert_values_equal_rel(r, std::sqrt(3.), 1.e-19);
    assert_values_equal_rel(theta, std::acos(-1. / std::sqrt(3.)), 1.e-19);
    assert_values_equal_rel(phi, -0.75 * M_PI, 1.e-19);

    const Direction direction(0., 0., 1.);
    const Point translated_point = point.translate(1., direction);
    assert_condition(translated_point.x() == point.x());
    assert_condition(translated_point.y() == point.y());
    assert_condition(translated_point.z() == point.z() + 1.);

    const Point rotated_point = point.rotate(direction, 0.5 * M_PI);
    ctm_warning("rotated point: %g %g %g", rotated_point.x(), rotated_point.y(),
                rotated_point.z());
    assert_values_equal_rel(rotated_point.x(), 1., 1.e-16);
    assert_values_equal_rel(rotated_point.y(), -1., 1.e-16);
    assert_condition(rotated_point.z() == point.z());
  }

  /// another special point
  {
    const Point point(0., 0., 0.);

    double r, theta, phi;
    point.spherical_coordinates(r, theta, phi);
    assert_condition(r == 0.);
    assert_condition(theta == 0.);
    assert_condition(phi == 0.);

    const Direction direction(0., 0., 1.);
    const Point translated_point = point.translate(1., direction);
    assert_condition(translated_point.x() == point.x());
    assert_condition(translated_point.y() == point.y());
    assert_condition(translated_point.z() == point.z() + 1.);

    const Point rotated_point = point.rotate(direction, 0.5 * M_PI);
    ctm_warning("rotated point: %g %g %g", rotated_point.x(), rotated_point.y(),
                rotated_point.z());
    assert_condition(rotated_point.x() == 0.);
    assert_condition(rotated_point.y() == 0.);
    assert_condition(rotated_point.z() == 0.);
  }

  /// a rotation similar to the one employed by random line functions
  {
    const Direction axis(0.5 * M_PI, M_PI / 6.);
    const Direction normal(M_PI / 3., -M_PI / 3.);
    RandomGenerator random_generator(42);

    for (uint_fast32_t i = 0; i < 1000000u; ++i) {
      const double r = std::sqrt(random_generator.get_uniform_random_double());
      const double phi =
          2. * M_PI * random_generator.get_uniform_random_double();
      const Point point(r * std::cos(phi), r * std::sin(phi), 0.);

      const Point rotated_point = point.rotate(axis, M_PI / 3.);

      const double nr = normal.nx() * rotated_point.x() +
                        normal.ny() * rotated_point.y() +
                        normal.nz() * rotated_point.z();
      assert_condition(std::abs(nr) < 1.e-15);
    }
  }

  return 0.;
}
