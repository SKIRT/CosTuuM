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
 * @file SizeBasedAlignmentDistribution.hpp
 *
 * @brief AlignmentDistribution based on the size of the particle.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SIZEBASEDALIGNMENTDISTRIBUTION_HPP
#define SIZEBASEDALIGNMENTDISTRIBUTION_HPP

#include "AlignmentDistribution.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "DisabledAlignmentOrientationDistribution.hpp"
#include "MishchenkoOrientationDistribution.hpp"

#include <cmath>

/**
 * @brief AlignmentDistribution based on the size of the particle.
 */
class SizeBasedAlignmentDistribution : public AlignmentDistribution {
private:
  /*! @brief Transition size above which dust grains are assumed to align with
   *  the magnetic field (in m). */
  const float_type _minimum_size;

  /*! @brief Orientation distribution for non-aligned grains. */
  OrientationDistribution *_non_aligned_orientation_distribution;

  /*! @brief Orientation distribution for oblate aligned grains. */
  OrientationDistribution *_oblate_aligned_orientation_distribution;

  /*! @brief Orientation distribution for prolate aligned grains. */
  OrientationDistribution *_prolate_aligned_orientation_distribution;

public:
  /**
   * @brief Constructor.
   *
   * @param minimum_size Transition size above which dust grains are assumed to
   * align with the magnetic field (in m).
   * @param aligned_orientation_distribution_type Type of orientation
   * distribution for aligned particles. Needs to be documented properly.
   * @param nmax Maximum order of spherical basis function expansions.
   * @param oblate_aligned_orientation_distribution Orientation distribution
   * for oblate aligned grains (if a custom type was chosen).
   * @param prolate_aligned_orientation_distribution Orientation distribution
   * for prolate aligned grains (if a custom type was chosen).
   */
  inline SizeBasedAlignmentDistribution(
      const float_type minimum_size,
      const int_fast32_t aligned_orientation_distribution_type,
      uint_fast32_t nmax,
      OrientationDistribution *oblate_aligned_orientation_distribution =
          nullptr,
      OrientationDistribution *prolate_aligned_orientation_distribution =
          nullptr)
      : _minimum_size(minimum_size),
        _non_aligned_orientation_distribution(nullptr),
        _oblate_aligned_orientation_distribution(nullptr),
        _prolate_aligned_orientation_distribution(nullptr) {

    _non_aligned_orientation_distribution =
        new OrientationDistribution(2 * nmax);
    _non_aligned_orientation_distribution->initialise();
    if (aligned_orientation_distribution_type == 0) {
      _oblate_aligned_orientation_distribution =
          new DavisGreensteinOrientationDistribution(2 * nmax, 2.);
      _prolate_aligned_orientation_distribution =
          new DavisGreensteinOrientationDistribution(2 * nmax, 0.5);
    } else if (aligned_orientation_distribution_type == 1) {
      // note that we hardcode the cos2beta values for now:
      //  - oblate: 3/5
      //  - prolate: 1/5
      _oblate_aligned_orientation_distribution =
          new MishchenkoOrientationDistribution(2 * nmax, 0.6);
      _prolate_aligned_orientation_distribution =
          new MishchenkoOrientationDistribution(2 * nmax, 0.2);
    } else if (aligned_orientation_distribution_type == 2) {
      _oblate_aligned_orientation_distribution =
          new DisabledAlignmentOrientationDistribution();
      _prolate_aligned_orientation_distribution =
          new DisabledAlignmentOrientationDistribution();
    } else if (aligned_orientation_distribution_type == 3) {
      ctm_assert(oblate_aligned_orientation_distribution != nullptr);
      _oblate_aligned_orientation_distribution =
          new OrientationDistribution(*oblate_aligned_orientation_distribution);
      ctm_assert(prolate_aligned_orientation_distribution != nullptr);
      _prolate_aligned_orientation_distribution = new OrientationDistribution(
          *prolate_aligned_orientation_distribution);
    } else {
      ctm_error("Unknown orientation distribution!");
    }
  }

  /**
   * @brief Destructor.
   */
  virtual ~SizeBasedAlignmentDistribution() {
    delete _non_aligned_orientation_distribution;
    delete _oblate_aligned_orientation_distribution;
    delete _prolate_aligned_orientation_distribution;
  }

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
                   const float_type axis_ratio) const {

    if (equal_volume_radius > _minimum_size) {
      if (axis_ratio > 1.) {
        return *_oblate_aligned_orientation_distribution;
      } else {
        return *_prolate_aligned_orientation_distribution;
      }
    } else {
      return *_non_aligned_orientation_distribution;
    }
  }
};

#endif // SIZEBASEDALIGNMENTDISTRIBUTION_HPP
