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

public:
  /**
   * @brief Constructor.
   *
   * @param compositions Composition parameter axis.
   * @param sizes Dust grain size parameter axis.
   * @param wavelengths Wavelength parameter axis.
   */
  inline ResultKey(const std::vector<int_fast32_t> &compositions,
                   const std::vector<float_type> &sizes,
                   const std::vector<float_type> &wavelengths)
      : _compositions(compositions), _sizes(sizes), _wavelengths(wavelengths) {}

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
   * @return Corresponding index in the result array.
   */
  inline size_t get_result_index(const size_t composition_index,
                                 const size_t size_index,
                                 const size_t wavelength_index) const {

    ctm_assert(composition_index < _compositions.size());
    ctm_assert(size_index < _sizes.size());
    ctm_assert(wavelength_index < _wavelengths.size());

    return composition_index * _sizes.size() * _wavelengths.size() +
           size_index * _wavelengths.size() + wavelength_index;
  }
};

#endif // RESULTKEY_HPP
