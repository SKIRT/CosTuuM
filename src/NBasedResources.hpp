/**
 * @file NBasedResources.hpp
 *
 * @brief Precomputed factors that only depend on the maximum order,
 * @f$n_{max}@f$.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef NBASEDRESOURCES_HPP
#define NBASEDRESOURCES_HPP

#include "Configuration.hpp"
#include "Matrix.hpp"
#include "QuickSchedWrapper.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Precomputed factors that only depend on the maximum order,
 * @f$n_{max}@f$.
 */
class NBasedResources : public Resource, public Task {
private:
  /*! @brief Precomputed factors @f$\sqrt{\frac{2n+1}{n(n+1)}}@f$ (array of size
   *  @f$n_{max}@f$). */
  std::vector<float_type> _dd;

  /*! @brief Precomputed factors
   *  @f$\frac{1}{2}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$
   *  (@f$n_{max}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _ann;

public:
  /**
   * @brief Constructor.
   *
   * @param maximum_order Maximum order of factors that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline NBasedResources(const uint_fast32_t maximum_order)
      : _dd(maximum_order, float_type(0.)), _ann(maximum_order, maximum_order) {
  }

  /**
   * @brief Compute the factors.
   */
  virtual void execute() {
    const uint_fast32_t nmax = _dd.size();
    for (uint_fast32_t ni = 0; ni < nmax; ++ni) {
      const float_type nn((ni + 2.) * (ni + 1.));
      const float_type d = sqrt((2. * (ni + 1.) + 1.) / nn);
      _dd[ni] = d;
      for (uint_fast32_t nj = 0; nj < ni + 1; ++nj) {
        const float_type ddd = 0.5 * d * _dd[nj];
        _ann(ni, nj) = ddd;
        _ann(nj, ni) = ddd;
      }
    }
  }

  /**
   * @brief Get the precomputed factor @f$\sqrt{\frac{2n+1}{n(n+1)}}@f$ for the
   * given value of @f$n@f$.
   *
   * @param n @f$n@f$.
   * @return @f$\sqrt{\frac{2n+1}{n(n+1)}}@f$.
   */
  inline float_type get_dd(const uint_fast32_t n) const {
    ctm_assert(n - 1 < _dd.size());
    // check that the resource was actually computed
    ctm_assert(_dd[0] != 0);
    return _dd[n - 1];
  }

  /**
   * @brief Get the precomputed factor
   * @f$\frac{1}{2}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$ for the given
   * value of @f$n@f$ and @f$n'@f$.
   *
   * @param n1 @f$n@f$.
   * @param n2 @f$n'@f$.
   * @return @f$\frac{1}{2}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$.
   */
  inline float_type get_ann(const uint_fast32_t n1,
                            const uint_fast32_t n2) const {
    ctm_assert(n1 - 1 < _dd.size());
    ctm_assert(n2 - 1 < _dd.size());
    // check that the resource was actually computed
    ctm_assert(_dd[0] != 0);
    return _ann(n1 - 1, n2 - 1);
  }
};

#endif // NBASEDRESOURCES_HPP
