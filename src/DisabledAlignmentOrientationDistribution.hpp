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
