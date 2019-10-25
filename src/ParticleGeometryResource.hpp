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
class ParticleGeometryResource : public Resource, public Task {
private:
  /*! @brief Equal volume sphere radius, @f$R_V@f$. */
  const float_type _R_V;

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
  const GaussBasedResources _quadrature_points;

public:
  /**
   * @brief Constructor.
   *
   * @param R_V Equal volume sphere radius, @f$R_V@f$.
   * @param axis_ratio Axis ratio of the spheroid, @f$d = \frac{a}{b}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param quadrature_points Gauss-Legendre quadrature points (to be computed
   * before this task is executed!).
   */
  inline ParticleGeometryResource(const float_type R_V,
                                  const float_type axis_ratio,
                                  const uint_fast32_t ngauss,
                                  const GaussBasedResources &quadrature_points)
      : _R_V(R_V), _axis_ratio(axis_ratio), _r(2 * ngauss, 0.),
        _r2(2 * ngauss, 0.), _dr_over_r(2 * ngauss, 0.),
        _quadrature_points(quadrature_points) {}

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
   * @brief Compute the factors.
   */
  virtual void execute() {
    SpecialFunctions::get_r_dr_spheroid(_quadrature_points.get_costhetas(),
                                        _R_V, _axis_ratio, _r2, _dr_over_r);
    for (uint_fast32_t ig = 0; ig < _r2.size(); ++ig) {
      _r[ig] = sqrt(_r2[ig]);
    }
  }

  /**
   * @brief Get the radius for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radius, @f$r@f$.
   */
  inline float_type get_r(const uint_fast32_t ig) const { return _r[ig]; }

  /**
   * @brief Get the radius squared for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radius squared, @f$r^2@f$.
   */
  inline float_type get_r2(const uint_fast32_t ig) const { return _r2[ig]; }

  /**
   * @brief Get the radial derivative for the given quadrature point.
   *
   * @param ig Index of the Gauss-Legendre quadrature point.
   * @return Corresponding radial derivative, @f$\frac{1}{r(\theta{})}
   * \frac{d}{d\theta{}} r(\theta{})@f$.
   */
  inline float_type get_dr_over_r(const uint_fast32_t ig) const {
    return _dr_over_r[ig];
  }
};

#endif // PARTICLEGEOMETRYRESOURCE_HPP
