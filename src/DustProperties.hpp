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
  /*! @brief Carbon grains. */
  DUSTGRAINTYPE_CARBON = 0,
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
