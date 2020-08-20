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
 * @file HEALPixRecorder.hpp
 *
 * @brief HEALPix unit sphere decomposition that can be used to bin directional
 * information.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEALPIXRECORDER_HPP
#define HEALPIXRECORDER_HPP

#include "Direction.hpp"

#include <cinttypes>
#include <fstream>
#include <vector>

/**
 * @brief HEALPix unit sphere decomposition that can be used to bin directional
 * information.
 */
template <uint_fast32_t _cell_size_> class HEALPixRecorder {
private:
  /*! @brief Order of the HEALPix grid. */
  const uint_fast8_t _order;

  /*! @brief Side length of the HEALPix grid. */
  const uint_fast32_t _Nside;

  /*! @brief Number of pixels in the polar region of the HEALPix grid. */
  const uint_fast32_t _Npixel_pole;

  /*! @brief Pixel values. */
  std::vector<double> _pixel_values;

public:
  /**
   * @brief Constructor.
   *
   * @param order Order of the HEALPix grid.
   */
  inline HEALPixRecorder(const uint_fast8_t order)
      : _order(order), _Nside(1u << _order),
        _Npixel_pole(2 * (_Nside * _Nside - _Nside)) {
    const uint_fast32_t Npixel = 12 * _Nside * _Nside;
    _pixel_values.resize(Npixel * _cell_size_, 0.);
  }

  /**
   * @brief Get the RING index of the HEALPix pixel that contains the given
   * direction.
   *
   * @param direction Direction.
   * @return Corresponding HEALPix pixel index in RING ordering.
   */
  inline uint_fast32_t ring_index(const Direction direction) const {

    double theta = direction.get_zenith_angle();
    double phi = direction.get_azimuth_angle();

    // the code below was mostly copied from healpix_base.cc in the official
    // HEALPix repository
    // it corresponds to the loc2pix function using the RING scheme
    const double z = std::cos(theta);
    const double za = std::abs(z);
    const double tt = std::fmod(2. * phi / M_PI, 4.);

    // below, the following conventions are used:
    //  i: pixel-in-ring index (horizontal pixel index)
    //  j: ring index (vertical pixel index)
    // first figure out if we are in the equatorial or polar region
    if (za <= 2. / 3.) {
      // equatorial region: all rings have 4*_Nside pixels
      const double t1 = _Nside * (0.5 + tt);
      const double t2 = 0.75 * _Nside * z;
      const int_fast32_t jp = static_cast<int_fast32_t>(std::floor(t1 - t2));
      const int_fast32_t jm = static_cast<int_fast32_t>(std::floor(t1 + t2));

      // we compute the ring index as the offset in the equatorial region
      // (j in [1, 2*_Nside+1])
      const uint_fast32_t j = _Nside + 1 + jp - jm;
      const bool kshift = 1 - (j & 1);
      const uint_fast32_t temp = jp + jm + kshift + 1 + 7 * _Nside;
      const uint_fast32_t i = (_order > 0) ? ((temp >> 1) & (4 * _Nside - 1))
                                           : ((temp >> 1) % (4 * _Nside));

      return _Npixel_pole + 4 * (j - 1) * _Nside + i;
    } else {
      // polar region: the number of pixels per ring depends on how far the ring
      // is from the pole
      const double tp = tt - static_cast<int>(tt);
      const double tmp =
          (za < 0.99) ? (_Nside * std::sqrt(3. * (1. - za)))
                      : (_Nside * std::sin(theta) / std::sqrt((1. + za) / 3.));

      const int_fast32_t jp = static_cast<int_fast32_t>(tp * tmp);
      const int_fast32_t jm = static_cast<int_fast32_t>((1. - tp) * tmp);

      // we compute the ring index as an offset from the pole (j in
      // [1, _Nside-1])
      const uint_fast32_t j = jp + jm + 1;
      const uint_fast32_t i = static_cast<uint_fast32_t>(tt * j);
      if (z < 0) {
        return 12 * _Nside * _Nside - 2 * j * (j + 1) + i;
      } else {
        return 2 * j * (j - 1) + i;
      }
    }
  }

  /**
   * @brief Bin the given data in the pixel that contains the given direction.
   *
   * @param direction Direction.
   * @param values Values to bin (array of length _cell_size_).
   */
  inline void bin(const Direction direction, const double *values) {
    const uint_fast32_t ibin = ring_index(direction);
    for (uint_fast32_t i = 0; i < _cell_size_; ++i) {
      _pixel_values[ibin * _cell_size_ + i] += values[i];
    }
  }

  /**
   * @brief Dump the recorded values to a binary file with the given name.
   *
   * @param filename Name of the binary file to create.
   */
  inline void dump(const std::string filename) const {
    std::ofstream file(filename);
    for (uint_fast32_t i = 0; i < _pixel_values.size(); ++i) {
      file.write(reinterpret_cast<const char *>(&_pixel_values[i]),
                 sizeof(double));
    }
    file.close();
  }
};

#endif // HEALPIXRECORDER_HPP
