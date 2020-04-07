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
 * @file GaussBasedResources.hpp
 *
 * @brief Precomputed factors that only depend on the number of Gauss-Legendre
 * quadrature points, @f$n_{GL}@f$.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GAUSSBASEDRESOURCES_HPP
#define GAUSSBASEDRESOURCES_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "QuickSchedWrapper.hpp"
#include "SpecialFunctions.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Precomputed factors that only depend on the number of Gauss-Legendre
 * quadrature points, @f$n_{GL}@f$.
 */
class GaussBasedResources : public Resource, public Task, public Computable {
private:
  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Precomputed factors @f$\cos(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _costheta;

  /*! @brief Precomputed factors @f$\frac{1}{\sin(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _sinthetainv;

  /*! @brief Precomputed factors @f$\frac{1}{\sin^2(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _sintheta2inv;

  /*! @brief Gauss-Legendre weights for the roots @f$\cos(\theta{})@f$ (array of
   *  size @f$2n_{GL}@f$). */
  std::vector<float_type> _weights;

public:
  /**
   * @brief Constructor.
   *
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline GaussBasedResources(const uint_fast32_t ngauss)
      : _ngauss(ngauss), _costheta(2 * ngauss, float_type(0.)),
        _sinthetainv(2 * ngauss, float_type(0.)),
        _sintheta2inv(2 * ngauss, float_type(0.)),
        _weights(2 * ngauss, float_type(0.)) {}

  virtual ~GaussBasedResources() {}

  /**
   * @brief Get the size in memory of a hypothetical GaussBasedResources object
   * with the given parameters.
   *
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t ngauss) {
    size_t size = sizeof(GaussBasedResources);
    // costheta
    size += 2 * ngauss * sizeof(float_type);
    // sinthetainv
    size += 2 * ngauss * sizeof(float_type);
    // sintheta2inv
    size += 2 * ngauss * sizeof(float_type);
    // weights
    size += 2 * ngauss * sizeof(float_type);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);
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
   * @return 0.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 0; }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    SpecialFunctions::get_gauss_legendre_points_and_weights(
        2 * _ngauss, _costheta, _weights);
    for (uint_fast32_t ig = 0; ig < _ngauss; ++ig) {
      const float_type this_sintheta2inv =
          1. / ((1. - _costheta[ig]) * (1. + _costheta[ig]));
      _sintheta2inv[ig] = this_sintheta2inv;
      _sintheta2inv[2 * _ngauss - ig - 1] = this_sintheta2inv;
      const float_type this_sinthetainv = sqrt(this_sintheta2inv);
      _sinthetainv[ig] = this_sinthetainv;
      _sinthetainv[2 * _ngauss - ig - 1] = this_sinthetainv;
    }

    ctm_assert_no_nans(_costheta, 2 * _ngauss);
    ctm_assert_no_nans(_weights, 2 * _ngauss);
    ctm_assert_no_nans(_sintheta2inv, 2 * _ngauss);
    ctm_assert_no_nans(_sinthetainv, 2 * _ngauss);

    make_available();
  }

  /**
   * @brief Get the computational cost of this task.
   *
   * @return Computational cost.
   */
  virtual int_fast32_t get_cost() const {
    return 208 * _ngauss * _ngauss + 5153 * _ngauss + 6035;
  }

  /**
   * @brief Get the precomputed @f$\cos(\theta{})@f$ value for the Gauss-
   * Legendre quadrature point with the given index.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Value of @f$\cos(\theta{})@f$ for that quadrature point.
   */
  inline float_type get_costheta(const uint_fast32_t ig) const {
    ctm_assert(ig < _costheta.size());
    // check that the resource was actually computed
    check_use();
    return _costheta[ig];
  }

  /**
   * @brief Get read-only access to the @f$\cos(\theta{})@f$ array.
   *
   * @return Reference to the @f$\cos(\theta{})@f$ array.
   */
  inline const std::vector<float_type> &get_costhetas() const {
    // check that the resource was actually computed
    check_use();
    return _costheta;
  }

  /**
   * @brief Get the precomputed @f$\frac{1}{\sin(\theta{})}@f$ value for the
   * Gauss-Legendre quadrature point with the given index.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Value of @f$\frac{1}{\sin(\theta{})}@f$ for that quadrature point.
   */
  inline float_type get_sinthetainv(const uint_fast32_t ig) const {
    ctm_assert(ig < _sinthetainv.size());
    // check that the resource was actually computed
    check_use();
    return _sinthetainv[ig];
  }

  /**
   * @brief Get the precomputed @f$\frac{1}{\sin^2(\theta{})}@f$ value for the
   * Gauss-Legendre quadrature point with the given index.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Value of @f$\frac{1}{\sin^2(\theta{})}@f$ for that quadrature
   * point.
   */
  inline float_type get_sintheta2inv(const uint_fast32_t ig) const {
    ctm_assert(ig < _sintheta2inv.size());
    // check that the resource was actually computed
    check_use();
    return _sintheta2inv[ig];
  }

  /**
   * @brief Get the precomputed weight for the Gauss-Legendre quadrature point
   * with the given index.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Weight for that quadrature point.
   */
  inline float_type get_weight(const uint_fast32_t ig) const {
    ctm_assert(ig < _weights.size());
    // check that the resource was actually computed
    check_use();
    return _weights[ig];
  }
};

#endif // GAUSSBASEDRESOURCES_HPP
