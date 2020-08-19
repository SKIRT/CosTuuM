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
 * @file RefractionEvent.hpp
 *
 * @brief Calculations for a refraction event at a surface between materials
 * with different refractive indices.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef REFRACTIONEVENT_HPP
#define REFRACTIONEVENT_HPP

#include "Error.hpp"

#include <cmath>

/**
 * @brief Possible polarisation states for a refracting ray.
 *
 * These are expressed w.r.t. the refraction plane, i.e. the plane through the
 * incident ray and the normal to the refraction surface.
 */
enum RefractionPolarisationState {
  /*! @brief Component perpendicular to the refraction plane. */
  S_POLARISED,
  /*! @brief Component parallel to the refraction plane. */
  P_POLARISED
};

/**
 * @brief Refraction event at a surface between materials with different
 * refractive indices.
 */
class RefractionEvent {
private:
  /*! @brief Reflection angle (in radians). */
  double _reflection_angle;

  /*! @brief Refraction angle (in radians). Set to a negative value if there is
   *  no refraction (total internal reflection). */
  double _refraction_angle;

  /*! @brief Reflection coefficients. */
  double _reflection_coefficient[2];

  /*! @brief Transmission coefficients. */
  double _transmission_coefficient[2];

public:
  /**
   * @brief Constructor.
   *
   * @param incident_angle Incident angle, @f$\theta{}_i@f$ (in radians).
   * @param n1 Refractive index of the material in which the incident ray is
   * travelling.
   * @param n2 Refractive index of the material on the other side of the
   * surface.
   */
  inline RefractionEvent(const double incident_angle, const double n1,
                         const double n2) {

    // sanity checks on the input values
    ctm_assert(incident_angle >= 0.);
    ctm_assert(incident_angle <= 0.5 * M_PI);
    ctm_assert(n1 > 0.);
    ctm_assert(n2 > 0.);

    // computing the reflection angle is trivial
    _reflection_angle = incident_angle;

    // compute the sine of the incident angle
    const double sintheta_i = std::sin(incident_angle);
    // use Snell's law to compute the sine of the refraction angle
    const double sintheta_t = n1 * sintheta_i / n2;

    // now check if we have refraction
    if (sintheta_t <= 1.) {
      _refraction_angle = std::asin(sintheta_t);
      const double costheta_i = std::cos(incident_angle);
      const double costheta_t = std::cos(_refraction_angle);
      // compute the reflection and transmission coefficients using
      // the Fresnel equations
      const double sqrtRs = (n1 * costheta_i - n2 * costheta_t) /
                            (n1 * costheta_i + n2 * costheta_t);
      const double sqrtRp = (n1 * costheta_t - n2 * costheta_i) /
                            (n1 * costheta_t + n2 * costheta_i);
      _reflection_coefficient[S_POLARISED] = sqrtRs * sqrtRs;
      _reflection_coefficient[P_POLARISED] = sqrtRp * sqrtRp;
      _transmission_coefficient[S_POLARISED] = 1. - _reflection_coefficient[0];
      _transmission_coefficient[P_POLARISED] = 1. - _reflection_coefficient[1];
    } else {
      // total internal reflection
      _refraction_angle = -1.;
      _reflection_coefficient[S_POLARISED] = 1.;
      _reflection_coefficient[P_POLARISED] = 1.;
      _transmission_coefficient[S_POLARISED] = 0.;
      _transmission_coefficient[P_POLARISED] = 0.;
    }
  }

  /**
   * @brief Get the reflection angle.
   *
   * @return Reflection angle (in radians).
   */
  inline double get_reflection_angle() const { return _reflection_angle; }

  /**
   * @brief Get the refraction angle.
   *
   * Note that this function returns a negative value if there is no refraction
   * (total internal reflection).
   *
   * @return Refraction angle (in radians).
   */
  inline double get_refraction_angle() const { return _refraction_angle; }

  /**
   * @brief Get the reflection coefficient for the given polarisation mode.
   *
   * @param polarisation_mode Valid RefractionPolarisationState.
   * @return Corresponding reflection coefficient.
   */
  inline double
  get_reflection_coefficient(const int_fast32_t polarisation_mode) const {
    return _reflection_coefficient[polarisation_mode];
  }

  /**
   * @brief Get the transmission coefficient for the given polarisation mode.
   *
   * @param polarisation_mode Valid RefractionPolarisationState.
   * @return Corresponding transmission coefficient.
   */
  inline double
  get_transmission_coefficient(const int_fast32_t polarisation_mode) const {
    return _transmission_coefficient[polarisation_mode];
  }
};

#endif // REFRACTIONEVENT_HPP
