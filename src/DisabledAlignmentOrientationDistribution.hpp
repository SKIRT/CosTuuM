/**
 * @file DisabledAlignmentOrientationDistribution.hpp
 *
 * @brief Dummy OrientationDistribution that does not compute any alignment
 * distribution, useful to expose the unaveraged T-matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DISABLEDALIGNMENTORIENTATIONDISTRIBUTION_HPP
#define DISABLEDALIGNMENTORIENTATIONDISTRIBUTION_HPP

#include "OrientationDistribution.hpp"

/**
 * @brief Dummy OrientationDistribution that does not compute any alignment
 * distribution, useful to expose the unaveraged T-matrix.
 */
class DisabledAlignmentOrientationDistribution
    : public OrientationDistribution {
public:
  /**
   * @brief Virtual function that contains an actual expression for the
   * orientation distribution function as a function of the zenith angle
   * @f$\beta{}@f$.
   *
   * The zenith angle @f$\beta{}@f$ represents the angle between the average
   * orientation of the rotation axis of the spheroidal dust particles and the
   * specific orientation for which we want to know the probability.
   * A distribution function @f$p(\beta{}) = \delta{}(\beta{})@f$ e.g. would
   * correspond to perfect alignment of all dust particles.
   *
   * The additional cosine and sine arguments are provided because they often
   * appear in more complicated expressions for @f$p(\beta{})@f$, and also
   * because they are computed in the integrand anyway.
   *
   * @param beta Zenith angle @f$\beta{} \in{} [0, \pi{}]@f$.
   * @param cosbeta Cosine of the zenith angle, @f$\cos(\beta{})@f$.
   * @param sinbeta Sine of the zenith angle, @f$\sin(\beta{})@f$.
   * @return Orientation distribution function for that value of
   * @f$\beta{}@f$, @f$p(\beta{})@f$.
   */
  virtual float_type operator()(const float_type beta, const float_type cosbeta,
                                const float_type sinbeta) const {

    ctm_error("Should not be used!");
  }

public:
  /**
   * @brief Constructor.
   */
  inline DisabledAlignmentOrientationDistribution()
      : OrientationDistribution(0) {}

  virtual ~DisabledAlignmentOrientationDistribution() {}

  /**
   * @brief Compute alignment for this distribution?
   *
   * @return True if the alignment should be computed.
   */
  virtual bool compute_alignment() const { return false; }
};

#endif // DISABLEDALIGNMENTORIENTATIONDISTRIBUTION_HPP
