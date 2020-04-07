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
 * @file DustProperties.hpp
 *
 * @brief General interface for classes containing dust properties.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DUSTPROPERTIES_HPP
#define DUSTPROPERTIES_HPP

#include "Configuration.hpp"

#include <complex>

/**
 * @brief Possible types of dust grains.
 */
enum DustGrainTypes {
  /*! @brief Carbon grains with their dipole moment oriented parallel to the
   *  incoming electromagnetic field. */
  DUSTGRAINTYPE_CARBON_PARALLEL = 0,
  /*! @brief Carbon grains with their dipole moment oriented perpendicular to
   *  the incoming electromagnetic field. */
  DUSTGRAINTYPE_CARBON_PERPENDICULAR,
  /*! @brief Silicon grains. */
  DUSTGRAINTYPE_SILICON,
  /*! @brief Number of dust grain types. */
  NUMBER_OF_DUSTGRAINTYPES
};

/**
 * @brief General interface for classes containing dust properties.
 */
class DustProperties {
public:
  virtual ~DustProperties() {}

  /**
   * @brief Get the refractive index for the given wavelength, grain size and
   * grain type.
   *
   * @param wavelength Wavelength of the incoming radiation (in m).
   * @param grain_size Size of the grain (in m).
   * @param grain_type Type of dust grain.
   * @return Complex refractive index of the grain for radiation at this
   * wavelength.
   */
  virtual std::complex<float_type>
  get_refractive_index(const float_type wavelength, const float_type grain_size,
                       const int_fast32_t grain_type) const = 0;
};

#endif // DUSTPROPERTIES_HPP
