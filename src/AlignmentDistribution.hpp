/**
 * @file AlignmentDistribution.hpp
 *
 * @brief Interface for alignment distributions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ALIGNMENTDISTRIBUTION_HPP
#define ALIGNMENTDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "OrientationDistribution.hpp"

#include <cmath>

/**
 * @brief Interface for alignment distributions.
 */
class AlignmentDistribution {
public:
  /**
   * @brief Virtual alignment distribution function.
   *
   * The default value assumes perfect alignment for grains with sizes larger
   * than @f$0.1@f$ @f$\mu{}@f$m and no alignment (random orientation) for
   * smaller grains.
   *
   * @param equal_volume_radius Size of the particle, parametrised as the radius
   * of a sphere with the same volume (in m).
   * @param axis_ratio Input axis ratio, @f$d = \frac{a}{b}@f$.
   * @param nmax Highest order of coefficient that is stored.
   * @return Pointer to an alignment distribution that can be used to compute
   * the alignment for this particle. Memory management for this pointer is
   * transferred to the caller.
   */
  virtual OrientationDistribution *
  operator()(const float_type equal_volume_radius, const float_type axis_ratio,
             const uint_fast32_t nmax) const {

    if (equal_volume_radius < 1.e-7) {
      return new OrientationDistribution(nmax);
    } else {
      return new DavisGreensteinOrientationDistribution(nmax, axis_ratio);
    }
  }
};

#endif // ALIGNMENTDISTRIBUTION_HPP
