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
 * @file testHEALPixRecorder.cpp
 *
 * @brief Unit test for the HEALPixRecorder class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "HEALPixRecorder.hpp"
#include "RandomGenerator.hpp"

#include <fstream>

/**
 * @brief Unit test for the HEALPixRecorder class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  HEALPixRecorder<1> healpix(4);
  // we write binary output, since healpy needs to use exactly the same angles
  // to produce exactly the same output
  std::ofstream ofile("test_healpix_recorder.dat");
  for (uint_fast32_t itheta = 0; itheta < 100u; ++itheta) {
    const double theta = (itheta + 0.5) * 0.01 * M_PI;
    for (uint_fast32_t iphi = 0; iphi < 200u; ++iphi) {
      const double phi = (iphi + 0.5) * 0.01 * M_PI;
      const Direction direction(theta, phi);
      const uint_fast32_t ring_index = healpix.ring_index(direction);
      ofile.write(reinterpret_cast<const char *>(&theta), sizeof(double));
      ofile.write(reinterpret_cast<const char *>(&phi), sizeof(double));
      ofile.write(reinterpret_cast<const char *>(&ring_index),
                  sizeof(uint_fast32_t));
    }
  }
  ofile.close();

  RandomGenerator random_generator;
  const double val[1] = {1.};
  for (uint_fast32_t i = 0; i < 100000u; ++i) {
    const Direction direction(
        std::acos(-1. + 2. * random_generator.get_uniform_random_double()),
        2. * M_PI * random_generator.get_uniform_random_double());
    healpix.bin(direction, val);
  }
  healpix.dump("test_healpix_output.dat");

  std::ofstream pfile("test_healpix_angles.dat");
  for (uint_fast32_t i = 0; i < healpix.get_number_of_pixels(); ++i) {
    const Direction direction = healpix.direction(i);
    const double theta = direction.get_zenith_angle();
    const double phi = direction.get_azimuth_angle();
    pfile.write(reinterpret_cast<const char *>(&theta), sizeof(double));
    pfile.write(reinterpret_cast<const char *>(&phi), sizeof(double));
  }
  pfile.close();

  return 0;
}
