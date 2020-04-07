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
 * @file ResultKey.hpp
 *
 * @brief Object that links input parameters to indices in an ordered results
 * array.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef RESULTKEY_HPP
#define RESULTKEY_HPP

#include "Configuration.hpp"
#include "Error.hpp"

#include <vector>

/**
 * @brief Interface for task-based computation results.
 */
class ResultKey {
private:
  /*! @brief Composition parameter axis. */
  const std::vector<int_fast32_t> _compositions;

  /*! @brief Dust grain size parameter axis. */
  const std::vector<float_type> _sizes;

  /*! @brief Wavelength parameter axis. */
  const std::vector<float_type> _wavelengths;

  /*! @brief Number of results per parameter set. */
  const uint_fast32_t _number_of_results;

public:
  /**
   * @brief Constructor.
   *
   * @param compositions Composition parameter axis.
   * @param sizes Dust grain size parameter axis.
   * @param wavelengths Wavelength parameter axis.
   * @param number_of_results Number of results per parameter set.
   */
  inline ResultKey(const std::vector<int_fast32_t> &compositions,
                   const std::vector<float_type> &sizes,
                   const std::vector<float_type> &wavelengths,
                   const uint_fast32_t number_of_results = 1)
      : _compositions(compositions), _sizes(sizes), _wavelengths(wavelengths),
        _number_of_results(number_of_results) {}

  /**
   * @brief Get the number of composition parameter values.
   *
   * @return Number of composition parameter values.
   */
  inline size_t composition_size() const { return _compositions.size(); }

  /**
   * @brief Get the number of dust grain size parameter values.
   *
   * @return Number of dust grain size parameter values.
   */
  inline size_t size_size() const { return _sizes.size(); }

  /**
   * @brief Get the number of wavelength parameter values.
   *
   * @return Number of wavelength parameter values.
   */
  inline size_t wavelength_size() const { return _wavelengths.size(); }

  /**
   * @brief Get the number of results per parameter set.
   *
   * @return Number of results per parameter set.
   */
  inline uint_fast32_t results_size() const { return _number_of_results; }

  /**
   * @brief Get the composition value corresponding to the given index.
   *
   * @param index Composition index.
   * @return Corresponding composition.
   */
  inline int_fast32_t get_composition(const size_t index) const {
    ctm_assert(index < _compositions.size());
    return _compositions[index];
  }

  /**
   * @brief Get the dust grain size value corresponding to the given index.
   *
   * @param index Dust grain size index.
   * @return Corresponding dust grain size (in m).
   */
  inline float_type get_size(const size_t index) const {
    ctm_assert(index < _sizes.size());
    return _sizes[index];
  }

  /**
   * @brief Get the wavelength value corresponding to the given index.
   *
   * @param index Wavelength index.
   * @return Corresponding wavelength (in m).
   */
  inline float_type get_wavelength(const size_t index) const {
    ctm_assert(index < _wavelengths.size());
    return _wavelengths[index];
  }

  /**
   * @brief Get the index in the result array corresponding to the given
   * parameter index values.
   *
   * The implementation needs to match the loop nesting used in
   * TaskManager::generate_tasks.
   *
   * @param composition_index Composition index.
   * @param size_index Dust grain size index.
   * @param wavelength_index Wavelength index.
   * @param result_index Index of a result in the list of all results for this
   * parameter set.
   * @return Corresponding index in the result array.
   */
  inline size_t get_result_index(const size_t composition_index,
                                 const size_t size_index,
                                 const size_t wavelength_index,
                                 const uint_fast32_t result_index = 0) const {

    ctm_assert(composition_index < _compositions.size());
    ctm_assert(size_index < _sizes.size());
    ctm_assert(wavelength_index < _wavelengths.size());
    ctm_assert(result_index < _number_of_results);

    return composition_index * _sizes.size() * _wavelengths.size() *
               _number_of_results +
           size_index * _wavelengths.size() * _number_of_results +
           wavelength_index * _number_of_results + result_index;
  }
};

#endif // RESULTKEY_HPP
