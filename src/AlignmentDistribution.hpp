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
#include "OrientationDistribution.hpp"

/**
 * @brief Interface for alignment distributions.
 */
class AlignmentDistribution {
public:
  virtual ~AlignmentDistribution() {}

  /**
   * @brief Get a reference to the alignment distribution function for a dust
   * grain with the given equal volume radius and axis ratio.
   *
   * @param equal_volume_radius Size of the particle, parametrised as the radius
   * of a sphere with the same volume (in m).
   * @param axis_ratio Input axis ratio, @f$d = \frac{a}{b}@f$.
   * @return Reference to the alignment distribution function for this dust
   * grain.
   */
  virtual const OrientationDistribution &
  get_distribution(const float_type equal_volume_radius,
                   const float_type axis_ratio) const = 0;
};

#endif // ALIGNMENTDISTRIBUTION_HPP
