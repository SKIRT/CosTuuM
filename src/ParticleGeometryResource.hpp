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
 * @file ParticleGeometryResource.hpp
 *
 * @brief Precomputed factors that depend on precomputed Gauss-Legendre factors,
 * and a specific choice of particle geometry.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PARTICLEGEOMETRYRESOURCE_HPP
#define PARTICLEGEOMETRYRESOURCE_HPP

#include "Configuration.hpp"
#include "GaussBasedResources.hpp"
#include "QuickSchedWrapper.hpp"
#include "SpecialFunctions.hpp"

#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif
/**
 * @brief Precomputed factors that depend on precomputed Gauss-Legendre factors,
 * and a specific choice of particle geometry.
 */
class ParticleGeometryResource : public Resource,
                                 public Task,
                                 public Computable {
private:
  /*! @brief Axis ratio of the spheroid, @f$d = \frac{a}{b}@f$. */
  const float_type _axis_ratio;

  /*! @brief Precomputed factors @f$r(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _r;

  /*! @brief Precomputed factors @f$r^2(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<float_type> _r2;

  /*! @brief Precomputed factors @f$\frac{1}{r(\theta{})}\frac{d}{d\theta{}}
   *  r(\theta{})@f$  (array of size @f$2n_{GL}@f$). */
  std::vector<float_type> _dr_over_r;

  /*! @brief Precomputed Gauss-Legendre quadrature points (read only). */
  const GaussBasedResources &_quadrature_points;

public:
  /**
   * @brief Constructor.
   *
   * @param axis_ratio Axis ratio of the spheroid, @f$d = \frac{a}{b}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param quadrature_points Gauss-Legendre quadrature points (to be computed
   * before this task is executed!).
   */
  inline ParticleGeometryResource(const float_type axis_ratio,
                                  const uint_fast32_t ngauss,
                                  const GaussBasedResources &quadrature_points)
      : _axis_ratio(axis_ratio), _r(2 * ngauss, 0.), _r2(2 * ngauss, 0.),
        _dr_over_r(2 * ngauss, 0.), _quadrature_points(quadrature_points) {}

  virtual ~ParticleGeometryResource() {}

  /**
   * @brief Get the size in memory of a hypothetical ParticleGeometryResource
   * object with the given parameters.
   *
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t ngauss) {
    size_t size = sizeof(ParticleGeometryResource);
    // r
    size += 2 * ngauss * sizeof(float_type);
    // r2
    size += 2 * ngauss * sizeof(float_type);
    // dr_over_r
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

    // read access
    quicksched.link_task_and_resource(*this, _quadrature_points, false);
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
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    SpecialFunctions::get_r_dr_spheroid(&_quadrature_points.get_costhetas()[0],
                                        _r2.size(), _axis_ratio, &_r2[0],
                                        &_dr_over_r[0]);
    for (uint_fast32_t ig = 0; ig < _r2.size(); ++ig) {
      _r[ig] = sqrt(_r2[ig]);
    }

    ctm_assert_no_nans(_r2, _r2.size());
    ctm_assert_no_nans(_r, _r.size());
    ctm_assert_no_nans(_dr_over_r, _dr_over_r.size());

    make_available();
  }

  /**
   * @brief Get the computational cost of this task.
   *
   * @return Computational cost.
   */
  virtual int_fast32_t get_cost() const { return 68 * _r2.size() + 2978; }

  /**
   * @brief Get the axis ratio for this geometry.
   *
   * @return Axis ratio, @f$d=\frac{a}{b}@f$.
   */
  inline float_type get_axis_ratio() const { return _axis_ratio; }

  /**
   * @brief Get the radius for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radius, @f$r@f$.
   */
  inline float_type get_r(const uint_fast32_t ig) const {
    ctm_assert(ig < _r.size());
    // check that the resource was actually computed
    check_use();
    return _r[ig];
  }

  /**
   * @brief Get the radius squared for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radius squared, @f$r^2@f$.
   */
  inline float_type get_r2(const uint_fast32_t ig) const {
    ctm_assert(ig < _r2.size());
    // check that the resource was actually computed
    check_use();
    return _r2[ig];
  }

  /**
   * @brief Get the radial derivative for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radial derivative, @f$\frac{1}{r(\theta{})}
   * \frac{d}{d\theta{}} r(\theta{})@f$.
   */
  inline float_type get_dr_over_r(const uint_fast32_t ig) const {
    ctm_assert(ig < _dr_over_r.size());
    // check that the resource was actually computed
    check_use();
    return _dr_over_r[ig];
  }
};

#endif // PARTICLEGEOMETRYRESOURCE_HPP
