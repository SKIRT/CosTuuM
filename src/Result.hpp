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
  RESULTTYPE_ABSORPTIONCOEFFICIENTS
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
};

#endif // RESULT_HPP
