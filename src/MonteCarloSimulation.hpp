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
 * @file MonteCarloSimulation.hpp
 *
 * @brief Monte Carlo simulation to determine the scattering properties of a
 * grain.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Grain.hpp"
#include "HEALPixRecorder.hpp"
#include "RefractionEvent.hpp"
#include "SphericalGrain.hpp"

#include <complex>

/**
 * @brief Monte Carlo simulation to determine the scattering properties of a
 * grain.
 */
class MonteCarloSimulation {
private:
  /*! @brief Grain. */
  Grain *_grain;

  /*! @brief Wavelength of the incoming radiation (in m). */
  double _wavelength;

  /*! @brief Equal volume radius of the grain (in m). */
  double _grain_size;

  /*! @brief Refractive index of the grain. */
  std::complex<double> _refractive_index;

  /*! @brief Direction of incident radiation. */
  Direction _direction;

  /*! @brief Number of integration points to use. */
  uint_fast64_t _number_of_points;

  /*! @brief Number of rays to shoot. */
  uint_fast64_t _number_of_rays;

  /*! @brief HEALPix recorder. */
  HEALPixRecorder<2u> _recorder;

  /*! @brief Random number generator. */
  RandomGenerator _random_generator;

public:
  /**
   * @brief Constructor.
   *
   * @param wavelength Wavelength of the incoming radiation (in m).
   * @param grain_size Equal volume radius of the grain (in m).
   * @param refractive_index Refractive index of the material.
   * @param direction Direction of incoming radiation.
   * @param number_of_points Number of integration points to use.
   * @param number_of_rays Number of rays to shoot.
   * @param healpix_order Order of the HEALPix grid used to record intensities.
   * @param random_seed Seed for the random number generator.
   */
  inline MonteCarloSimulation(const double wavelength, const double grain_size,
                              const std::complex<double> refractive_index,
                              const Direction direction,
                              const uint_fast64_t number_of_points,
                              const uint_fast64_t number_of_rays,
                              const uint_fast8_t healpix_order,
                              const int_fast32_t random_seed = 42)
      : _grain(new SphericalGrain()), _wavelength(wavelength),
        _grain_size(grain_size), _refractive_index(refractive_index),
        _direction(direction), _number_of_points(number_of_points),
        _number_of_rays(number_of_rays), _recorder(healpix_order),
        _random_generator(random_seed) {}

  /**
   * @brief Destructor.
   */
  inline ~MonteCarloSimulation() { delete _grain; }

  /**
   * @brief Run the Monte Carlo simulation.
   */
  inline void run() {

    const double k = 2. * M_PI / _wavelength;
    const double norm =
        0.25 * k * k * _grain_size * _grain_size * _grain_size * _grain_size;
    // add the central part of the diffraction pattern
    double values[2] = {0., norm};
    _recorder.bin(_direction, values);
    //    for (uint_fast32_t ipix = 0; ipix < _recorder.get_number_of_pixels();
    //         ++ipix) {
    //      const Direction direction = _recorder.direction(ipix);
    //      std::complex<double> integral(0., 0.);
    //      for (uint_fast32_t i = 0; i < _number_of_points; ++i) {
    //        const Point random_point =
    //        _grain->generate_random_cross_section_point(
    //            _direction, _random_generator);
    //        const Point point = random_point.scale(_grain_size);
    //        const double exponent =
    //            k * std::sin(direction.get_zenith_angle()) *
    //            (point.x() * std::cos(direction.get_azimuth_angle()) +
    //             point.y() * std::sin(direction.get_azimuth_angle()));
    //        const std::complex<double> integrand(std::cos(exponent),
    //                                             -std::sin(exponent));
    //        integral += integrand;
    //      }
    //      integral /= _number_of_points;
    //      values[1] = norm*std::norm(integral);
    //      _recorder.bin(direction, values);
    //    }

    for (uint_fast64_t iray = 0; iray < _number_of_rays; ++iray) {
      double I[2] = {M_PI * _grain_size * _grain_size / _number_of_rays, 0.};
      const Line line =
          _grain->generate_random_line(_direction, _random_generator);
      const IntersectionEvent intersection_event =
          _grain->get_intersection(line, false);
      const Point intersection_point =
          intersection_event.get_intersection_point();
      const Direction intersection_direction(intersection_point.x(),
                                             intersection_point.y(),
                                             intersection_point.z());
      _recorder.bin(intersection_direction, I);
      I[0] = 0.;
      // we need to divide the intensities by the recorder pixel size to
      // normalise them properly
      I[1] = M_PI * _grain_size * _grain_size /
             (_number_of_rays * 4. * M_PI / _recorder.get_number_of_pixels());
      const Direction normal = intersection_event.get_normal();
      const Direction reverse_direction = _direction.reverse();
      const double incident_angle = reverse_direction.angle(normal);
      const RefractionEvent refraction_event(incident_angle, 1.,
                                             _refractive_index.real());
      const double reflection_angle = refraction_event.get_reflection_angle();
      const Direction scatter_plane_normal = normal.cross(_direction);
      const Direction out = reverse_direction.rotate_perpendicular(
          scatter_plane_normal, incident_angle + reflection_angle);
      I[1] *= refraction_event.get_average_reflection_coefficient();
      _recorder.bin(out, I);
    }
  }

  /**
   * @brief Output the result of the Monte Carlo simulation.
   *
   * @param name Name of the binary file to dump.
   */
  inline void output(const std::string name) const { _recorder.dump(name); }
};

#endif // MONTECARLOSIMULATION_HPP
