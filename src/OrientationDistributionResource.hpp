/**
 * @file OrientationDistributionResource.hpp
 *
 * @brief Resource that wraps around an orientation distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ORIENTATIONDISTRIBUTIONRESOURCE_HPP
#define ORIENTATIONDISTRIBUTIONRESOURCE_HPP

#include "OrientationDistribution.hpp"
#include "QuickSchedWrapper.hpp"

/**
 * @brief Resource that wraps around an orientation distribution.
 */
class OrientationDistributionResource : public Resource {
private:
  /*! @brief Pointer to the wrapped orientation distribution. */
  OrientationDistribution *_orientation_distribution;

public:
  /**
   * @brief Constructor.
   *
   * @param orientation_distribution Pointer to the wrapped orientation
   * distribution.
   */
  inline OrientationDistributionResource(
      OrientationDistribution *orientation_distribution)
      : _orientation_distribution(orientation_distribution) {}

  virtual ~OrientationDistributionResource() {}

  /**
   * @brief Get the maximum available coefficient in the distribution.
   *
   * @return Maximum order of available coefficients.
   */
  inline uint_fast32_t get_maximum_order() const {
    return _orientation_distribution->get_maximum_order();
  }

  /**
   * @brief Get the coefficient for the given order.
   *
   * @param n Order @f$n@f$.
   * @return Corresponding coefficient.
   */
  inline float_type get_coefficient(const uint_fast32_t n) const {
    return _orientation_distribution->get_coefficient(n);
  }
};

#endif // ORIENTATIONDISTRIBUTIONRESOURCE_HPP
