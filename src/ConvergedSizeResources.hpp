/**
 * @file ConvergedSizeResources.hpp
 *
 * @brief Maximum order and number of Gauss-Legendre quadrature points for a
 * converged T-matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CONVERGEDSIZERESOURCES_HPP
#define CONVERGEDSIZERESOURCES_HPP

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
 * @brief Maximum order and number of Gauss-Legendre quadrature points for a
 * converged T-matrix.
 */
class ConvergedSizeResources : public Resource {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  uint_fast32_t _ngauss;

public:
  /**
   * @brief Empty constructor.
   */
  inline ConvergedSizeResources() : _nmax(0), _ngauss(0) {}

  /**
   * @brief Get the maximum order.
   *
   * @return Maximum order, @f$n_{max}@f$.
   */
  inline uint_fast32_t get_nmax() const { return _nmax; }

  /**
   * @brief Get the number of Gauss-Legendre quadrature points.
   *
   * @return Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline uint_fast32_t get_ngauss() const { return _ngauss; }
};

#endif // CONVERGEDSIZERESOURCES_HPP
