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
 * @file SphericalGrain.hpp
 *
 * @brief Grain implementation for spherical grains.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPHERICALGRAIN_HPP
#define SPHERICALGRAIN_HPP

#include "Grain.hpp"
#include "SphericalGrainProjection.hpp"

/**
 * @brief Grain implementation for spherical grains.
 */
class SphericalGrain : public Grain {
public:
  /**
   * @brief Get the projection of this grain perpendicular to the given
   * direction.
   *
   * @param direction Direction.
   * @return Projection perpendicular to this direction.
   */
  virtual GrainProjection get_projection(const Direction direction) const {
    return SphericalGrainProjection();
  }
};

#endif // SPHERICALGRAIN_HPP
