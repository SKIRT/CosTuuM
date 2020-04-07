/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
