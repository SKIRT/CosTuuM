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
 * @file testRefractionEvent.cpp
 *
 * @brief Unit test for the RefractionEvent class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "RefractionEvent.hpp"

/**
 * @brief Unit test for the RefractionEvent class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  /// normal refraction event
  {
    const RefractionEvent event(M_PI / 6., 1., 2.);

    ctm_warning("theta_r: %g", event.get_reflection_angle());
    ctm_warning("theta_t: %g", event.get_refraction_angle());
    ctm_warning("R: %g %g", event.get_reflection_coefficient(S_POLARISED),
                event.get_reflection_coefficient(P_POLARISED));
    ctm_warning("T: %g %g", event.get_transmission_coefficient(S_POLARISED),
                event.get_transmission_coefficient(P_POLARISED));

    // sin(theta_i) = sin(30 deg) = 1/2
    // sin(theta_t) = 1/2 * sin(theta_i) = 1/4
    // cos(theta_t) = sqrt(1 - sin(theta_t)) = sqrt(15)/4
    // cos(theta_i) = cos(30 deg) = sqrt(3)/2
    // hence:
    const double Rs =
        (std::sqrt(3.) - std::sqrt(15.)) / (std::sqrt(3.) + std::sqrt(15.)) *
        (std::sqrt(3.) - std::sqrt(15.)) / (std::sqrt(3.) + std::sqrt(15.));
    const double Rp = (std::sqrt(15.) - 4. * std::sqrt(3.)) /
                      (std::sqrt(15.) + 4. * std::sqrt(3.)) *
                      (std::sqrt(15.) - 4. * std::sqrt(3.)) /
                      (std::sqrt(15.) + 4. * std::sqrt(3.));

    assert_condition(event.get_reflection_angle() == M_PI / 6.);
    assert_condition(event.get_refraction_angle() == std::asin(0.25));
    assert_values_equal_rel(event.get_reflection_coefficient(S_POLARISED), Rs,
                            1.e-15);
    assert_values_equal_rel(event.get_reflection_coefficient(P_POLARISED), Rp,
                            1.e-15);
    assert_values_equal_rel(event.get_transmission_coefficient(S_POLARISED),
                            1. - Rs, 1.e-15);
    assert_values_equal_rel(event.get_transmission_coefficient(P_POLARISED),
                            1. - Rp, 1.e-15);
  }

  /// total internal reflection
  {
    const RefractionEvent event(M_PI / 3., 2., 1.);

    ctm_warning("theta_r: %g", event.get_reflection_angle());
    ctm_warning("theta_t: %g", event.get_refraction_angle());
    ctm_warning("R: %g %g", event.get_reflection_coefficient(S_POLARISED),
                event.get_reflection_coefficient(P_POLARISED));
    ctm_warning("T: %g %g", event.get_transmission_coefficient(S_POLARISED),
                event.get_transmission_coefficient(P_POLARISED));

    assert_condition(event.get_reflection_angle() == M_PI / 3.);
    assert_condition(event.get_refraction_angle() < 0.);
    assert_condition(event.get_reflection_coefficient(S_POLARISED) == 1.);
    assert_condition(event.get_reflection_coefficient(P_POLARISED) == 1.);
    assert_condition(event.get_transmission_coefficient(S_POLARISED) == 0.);
    assert_condition(event.get_transmission_coefficient(P_POLARISED) == 0.);
  }

  return 0;
}
