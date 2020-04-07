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
 * @file AlignmentAverageTask.hpp
 *
 * @brief Task that averages a T-matrix over a given alignment distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ALIGNMENTAVERAGETASK_HPP
#define ALIGNMENTAVERAGETASK_HPP

#include "OrientationDistribution.hpp"
#include "QuickSchedWrapper.hpp"
#include "TMatrixResource.hpp"

/**
 * @brief Task that averages a T-matrix over a given alignment distribution.
 */
class AlignmentAverageTask : public Task {
private:
  /*! @brief Orientation distribution. */
  const OrientationDistribution &_orientation_distribution;

  /*! @brief Input T-matrix. */
  const TMatrixResource &_input_Tmatrix;

  /*! @brief Auxiliary space manager used to obtain space to store intermediate
   *  calculations. */
  TMatrixAuxiliarySpaceManager &_aux_manager;

  /*! @brief Output T-matrix. */
  TMatrixResource &_output_Tmatrix;

public:
  /**
   * @brief Constructor.
   *
   * @param orientation_distribution Orientation distribution.
   * @param input_Tmatrix Input T-matrix.
   * @param aux_manager Auxiliary space manager used to obtain space to store
   * intermediate calculations.
   * @param output_Tmatrix Output T-matrix.
   */
  inline AlignmentAverageTask(
      const OrientationDistribution &orientation_distribution,
      const TMatrixResource &input_Tmatrix,
      TMatrixAuxiliarySpaceManager &aux_manager,
      TMatrixResource &output_Tmatrix)
      : _orientation_distribution(orientation_distribution),
        _input_Tmatrix(input_Tmatrix), _aux_manager(aux_manager),
        _output_Tmatrix(output_Tmatrix) {}

  virtual ~AlignmentAverageTask() {}

  /**
   * @brief Get the size in memory of a hypothetical AlignmentAverageTask object
   * with the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(AlignmentAverageTask);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _input_Tmatrix, true);

    // read access
    quicksched.link_task_and_resource(*this, _output_Tmatrix, false);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 1; }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t nmax = _input_Tmatrix.get_nmax();

    // first make sure the _nmax value for the new T-matrix is set
    _output_Tmatrix._nmax = nmax;

    // get the array in which to temporarily store Clebsch-Gordan coefficients
    std::vector<float_type> &CGcoeff =
        _aux_manager.get_space(thread_id)._clebsch_gordan;

    if (_orientation_distribution.compute_alignment()) {
      // check that we precomputed enough expansion coefficients in the
      // expansion of the orientation distribution
      ctm_assert(_orientation_distribution.get_maximum_order() >= 2 * nmax);

      // loop over all T matrix elements:
      //  - outer loop over m
      //  - inner loops over n1 and n2
      for (uint_fast32_t m = 0; m < nmax + 1; ++m) {
        // msign = (-1)^m
        const int_fast8_t msign = (m % 2 == 0) ? 1 : -1;
        // the minimum n value depends on the value of m
        const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));
        // inner loops
        for (uint_fast32_t n1 = nmin; n1 < nmax + 1; ++n1) {
          for (uint_fast32_t n2 = nmin; n2 < nmax + 1; ++n2) {
            // compute a new value for the T matrix coefficients T_{mn_1mn_2}
            std::complex<float_type> Tn1n2[2][2];
            // std::abs does not seem to work well in this context when quad
            // precision is activated, so we just mimic it ourselves
            const uint_fast32_t n12min = (n1 > n2) ? n1 - n2 : n2 - n1;
            const uint_fast32_t n12max = n1 + n2;
            const uint_fast32_t M = std::min(n1, n2);
            // first summation in Mishchenko (1991), eq. 3.27
            for (uint_fast32_t N = n12min; N < n12max + 1; ++N) {
              // now that we have n1, n2 and N, we can compute the
              // Clebsch-Gordan coefficients for this element of the sum
              //              const std::vector<float_type> CGcoeff =
              //                  SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(
              //                      n1, n2, N);
              SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(
                  n1, n2, N, CGcoeff);
              // we can also get the expansion coefficient
              const float_type pN =
                  _orientation_distribution.get_coefficient(N);
              // second summation in Mishchenko (1991), eq. 3.27
              for (uint_fast32_t m1 = 0; m1 < M + 1; ++m1) {
                // all elements m1=/=0 contribute twice because of symmetry
                float_type factor(2.);
                if (m1 == 0) {
                  factor = 1.;
                }
                // m1sign = (-1)^m1
                const int_fast8_t m1sign = (m1 % 2 == 0) ? 1 : -1;
                const int_fast8_t mm1sign = m1sign * msign;
                ctm_assert(M + m < CGcoeff.size());
                ctm_assert(M + m1 < CGcoeff.size());
                // prefactor for all terms in this part of the sum
                const float_type CGfac =
                    mm1sign * factor * pN * CGcoeff[M + m] * CGcoeff[M + m1];
                // depending on the sign of N + n1 + n2, odd/even terms in the
                // sum cancel out because of symmetry
                if ((N + n12max) % 2 == 0) {
                  Tn1n2[0][0] += CGfac * _input_Tmatrix(0, n1, m1, 0, n2, m1);
                  Tn1n2[1][1] += CGfac * _input_Tmatrix(1, n1, m1, 1, n2, m1);
                } else {
                  Tn1n2[0][1] += CGfac * _input_Tmatrix(0, n1, m1, 1, n2, m1);
                  Tn1n2[1][0] += CGfac * _input_Tmatrix(1, n1, m1, 0, n2, m1);
                }
              }
            }
            // now set the new values for T_{mn_1mn_2}
            _output_Tmatrix.set_element(0, n1, m, 0, n2, m, Tn1n2[0][0]);
            _output_Tmatrix.set_element(0, n1, m, 1, n2, m, Tn1n2[0][1]);
            _output_Tmatrix.set_element(1, n1, m, 0, n2, m, Tn1n2[1][0]);
            _output_Tmatrix.set_element(1, n1, m, 1, n2, m, Tn1n2[1][1]);
          }
        }
      }
    } else {
      // simply copy the T-matrix elements from the input to the output matrix
      for (uint_fast32_t m = 0; m < nmax + 1; ++m) {
        // the minimum n value depends on the value of m
        const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));
        // inner loops
        for (uint_fast32_t n1 = nmin; n1 < nmax + 1; ++n1) {
          for (uint_fast32_t n2 = nmin; n2 < nmax + 1; ++n2) {
            _output_Tmatrix.set_element(0, n1, m, 0, n2, m,
                                        _input_Tmatrix(0, n1, m, 0, n2, m));
            _output_Tmatrix.set_element(0, n1, m, 1, n2, m,
                                        _input_Tmatrix(0, n1, m, 1, n2, m));
            _output_Tmatrix.set_element(1, n1, m, 0, n2, m,
                                        _input_Tmatrix(1, n1, m, 0, n2, m));
            _output_Tmatrix.set_element(1, n1, m, 1, n2, m,
                                        _input_Tmatrix(1, n1, m, 1, n2, m));
          }
        }
      }
    }
  }
};

#endif // ALIGNMENTAVERAGETASK_HPP
