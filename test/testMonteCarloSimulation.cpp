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
 * @file testMonteCarloSimulation.cpp
 *
 * @brief Unit test for the MonteCarloSimulation class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "MonteCarloSimulation.hpp"

/**
 * @brief Unit test for the MonteCarloSimulation class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  MonteCarloSimulation simulation(std::complex<double>(2., 0.),
                                  Direction(0.5 * M_PI, 0.), 1000000u, 4, 42);
  simulation.run();
  simulation.output("test_monte_carlo_simulation.dat");

  return 0;
}
