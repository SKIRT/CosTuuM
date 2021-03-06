/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Assert.hpp
 *
 * @brief Custom assert macros
 *
 * This file was originally part of the photoionization code CMacIonize
 * (https://github.com/bwvdnbro/CMacIonize).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ASSERT_HPP
#define ASSERT_HPP

#include "Error.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

/*! @brief Assert that the given condition is true and throw an error when it is
 *  not. */
#define assert_condition(condition)                                            \
  if (!(condition)) {                                                          \
    ctm_error("Assertion failed (%s)!", #condition);                           \
  }

/*! @brief Assert that the given condition is true and throw an error with the
 *  given additional message if it is not. */
#define assert_condition_message(condition, message, ...)                      \
  {                                                                            \
    if (!(condition)) {                                                        \
      ctm_error("Assertion failed (%s): " message "!", #condition,             \
                ##__VA_ARGS__);                                                \
    }                                                                          \
  }

/*! @brief Assert that the given values are equal up to the given relative and
 *  absolute tolerance level. This is useful for comparing floating point values
 *  with round off error. */
#define assert_values_equal_tol(a, b, tol)                                     \
  if (std::abs((a) - (b)) > tol &&                                             \
      std::abs((a) - (b)) > tol * std::abs((a) + (b))) {                       \
    ctm_error("Assertion failed: %s (%g) != %s (%g)", #a, a, #b, b);           \
  }

/*! @brief Assert that the given values are equal up to the given relative
 *  tolerance level. This is useful for comparing floating point values with
 *  round off error. */
#define assert_values_equal_rel(a, b, tol)                                     \
  if (std::abs((a) - (b)) > tol * std::abs((a) + (b))) {                       \
    ctm_error(                                                                 \
        "Assertion failed: %s (%g) != %s (%g) (relative_difference: %g)", #a,  \
        a, #b, b, std::abs((a) - (b)) / std::abs((a) + (b)));                  \
  }

/*! @brief Assert that the given values are equal up to a pre-defined relative
 *  and absolute tolerance level of 1.e-4. This is useful for comparing floating
 *  point values with round off error. */
#define assert_values_equal(a, b) assert_values_equal_tol(a, b, 1.e-4)

#endif // ASSERT_HPP
