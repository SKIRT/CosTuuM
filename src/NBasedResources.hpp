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
class NBasedResources : public Resource, public Task, public Computable {
private:
  /*! @brief Precomputed factors
   *  @f$\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$
   *  (@f$n_{max}\times{}n_{max}@f$ matrix). */
  Matrix<float_type> _ann;

  /*! @brief Precomputed factors
   *  @f$i^{n'-n-1}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$
   *  (@f$n_{max}\times{}n_{max}@f$ matrix). */
  Matrix<std::complex<float_type>> _cnn;

public:
  /**
   * @brief Constructor.
   *
   * @param maximum_order Maximum order of factors that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline NBasedResources(const uint_fast32_t maximum_order)
      : _ann(maximum_order, maximum_order), _cnn(maximum_order, maximum_order) {
  }

  /**
   * @brief Get the size in memory of a hypothetical NBasedResources object with
   * the given parameters.
   *
   * @param maximum_order Maximum order of the hypothetical object.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t maximum_order) {
    size_t size = sizeof(NBasedResources);
    // ann
    size += maximum_order * maximum_order;
    // cnn
    size += 2 * maximum_order * maximum_order;
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
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    const uint_fast32_t nmax = _ann.get_number_of_rows();
    const std::complex<float_type> icompl(0., 1.);
    std::complex<float_type> icomp_pow_nn = icompl;
    for (uint_fast32_t ni = 0; ni < nmax; ++ni) {
      const float_type nni((ni + 2.) * (ni + 1.));
      const float_type sqrtni = sqrt((2. * (ni + 1.) + 1.) / nni);
      ctm_assert_not_nan(sqrtni);
      std::complex<float_type> icomp_pow_m_n_m_1(-1.);
      for (uint_fast32_t nj = 0; nj < nmax; ++nj) {
        const float_type nnj((nj + 2.) * (nj + 1.));
        const float_type sqrtnj = sqrt((2. * (nj + 1.) + 1.) / nnj);
        ctm_assert_not_nan(sqrtnj);
        _ann(ni, nj) = sqrtni * sqrtnj;
        ctm_assert_not_nan(_ann(ni, nj));
        _cnn(ni, nj) = icomp_pow_m_n_m_1 * icomp_pow_nn * _ann(ni, nj);
        icomp_pow_m_n_m_1 /= icompl;
      }
      icomp_pow_nn *= icompl;
    }
    make_available();
  }

  /**
   * @brief Get the precomputed factor
   * @f$\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$ for the given
   * value of @f$n@f$ and @f$n'@f$.
   *
   * @param n1 @f$n@f$.
   * @param n2 @f$n'@f$.
   * @return @f$\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$.
   */
  inline float_type get_ann(const uint_fast32_t n1,
                            const uint_fast32_t n2) const {
    ctm_assert(n1 - 1 < _ann.get_number_of_rows());
    ctm_assert(n2 - 1 < _ann.get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _ann(n1 - 1, n2 - 1);
  }

  /**
   * @brief Get the precomputed factor
   * @f$i^{n'-n-1}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$ for the given
   * value of @f$n@f$ and @f$n'@f$.
   *
   * @param n1 @f$n@f$.
   * @param n2 @f$n'@f$.
   * @return @f$i^{n'-n-1}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$.
   */
  inline std::complex<float_type> get_cnn(const uint_fast32_t n1,
                                          const uint_fast32_t n2) const {
    ctm_assert(n1 - 1 < _cnn.get_number_of_rows());
    ctm_assert(n2 - 1 < _cnn.get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _cnn(n1 - 1, n2 - 1);
  }
};

#endif // NBASEDRESOURCES_HPP
