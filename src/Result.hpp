/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Result.hpp
 *
 * @brief Interface for task-based computation results.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RESULT_HPP
#define RESULT_HPP

#include "Configuration.hpp"
#include "QuickSchedWrapper.hpp"

/**
 * @brief Types of results.
 */
enum ResultType {
  /*! @brief Absorption coefficient result. */
  RESULTTYPE_ABSORPTIONCOEFFICIENTS = 0,
  /*! @brief Extinction coefficient result. */
  RESULTTYPE_EXTINCTIONCOEFFICIENTS
};

/**
 * @brief Interface for task-based computation results.
 */
class Result : public Resource {
private:
  /*! @brief Composition parameter value for the result. */
  const int_fast32_t _composition;

  /*! @brief Particle size parameter value for the result. */
  const float_type _size;

  /*! @brief Wavelength value for the result. */
  const float_type _wavelength;

  /*! @brief Type of result. */
  const int_fast32_t _type;

public:
  /**
   * @brief Constructor.
   *
   * @param composition Composition parameter value for the result.
   * @param size Particle size parameter value for the result (in m).
   * @param wavelength Wavelength value for the result (in m).
   * @param type Type of result.
   */
  inline Result(const int_fast32_t composition, const float_type size,
                const float_type wavelength, const int_fast32_t type)
      : _composition(composition), _size(size), _wavelength(wavelength),
        _type(type) {

    (void)_size;
    (void)_composition;
    (void)_wavelength;
    (void)_type;
  }

  /**
   * @brief Get the composition for which this result was computed.
   *
   * @return Composition.
   */
  inline int_fast32_t get_composition() const { return _composition; }

  /**
   * @brief Get the size for which this result was computed.
   *
   * @return Size (in m).
   */
  inline float_type get_size() const { return _size; }

  /**
   * @brief Get the wavelength for which this result was computed.
   *
   * @return Wavelength (in m).
   */
  inline float_type get_wavelength() const { return _wavelength; }

  /**
   * @brief Get the type of result stored in this object.
   *
   * @return Result type.
   */
  inline int_fast32_t get_type() const { return _type; }
};

#endif // RESULT_HPP
