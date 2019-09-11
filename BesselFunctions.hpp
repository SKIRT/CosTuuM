/**
 * @file BesselFunctions.hpp
 *
 * @brief Class namespace that contains template functions to compute general
 * spherical Bessel functions of the first and second kind and a variant of
 * their derivatives, for real and complex input values.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include <cinttypes>
#include <cmath>
#include <complex>

/*! @brief Starting order for the backwards recurrence algorithm for the
 *  spherical Bessel functions of the first kind. The higher this number, the
 *  higher the accuracy of the Bessel functions of the first kind, but also
 *  the slower the function becomes. */
#define BESSELFUNCTIONS_NMAX 800u

/**
 * @brief Class namespace that contains template functions to compute general
 * spherical Bessel functions of the first and second kind and a variant of
 * their derivatives, for real and complex input values.
 */
class BesselFunctions {

public:
  /**
   * @brief Spherical Bessel function of the second kind for general real or
   * complex numbers that returns all spherical Bessel functions and a special
   * form of their first derivatives up to the given maximum order.
   *
   * Since this method is not part of the standard C++ library, we have to
   * implement it ourselves. We make use of the following two recursion
   * relations for spherical Bessel functions:
   * @f[
   *    y_{n-1}(z) + y_{n+1}(z) = \frac{2n + 1}{z} y_n(z),
   * @f]
   * and
   * @f[
   *    \frac{d}{dz} y_n(z) = y_{n-1}(z) - \frac{n+1}{z} y_n(z),
   * @f]
   * combined with our knowledge of the first and second order spherical Bessel
   * function of the second kind:
   * @f[
   *    y_1(z) = -\frac{\cos(z)}{z^2} - \frac{\sin(z)}{z},
   * @f]
   * and
   * @f[
   *    y_2(z) = \left(-\frac{3}{z^2} + 1\right) \frac{\cos(z)}{z}
   *             - 3 \frac{\sin(z)}{z^2}.
   * @f]
   *
   * We use a forward recursion algorithm to determine all values of
   * @f$y_n(z)@f$ and its derivative.
   *
   * Since we do not require the actual derivative, but rather the expression
   * @f[
   *    \frac{(zy_n(z))'}{z} = y_n'(z) + \frac{y_n(z)}{z},
   * @f]
   * we make use of the fact that this expression can be computed very easily
   * during the recursion step and skip the explicit computation of the actual
   * derivative.
   *
   * @param nmax Maximum order to compute (we compute all @f$y_n(z)@f$ for
   * @f$n\in{}[1,n_{max}]@f$).
   * @param z Input values.
   * @param y Array to store the Bessel function values in (of size nmax).
   * @param dy Array to store the first derivatives in (of size nmax).
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void spherical_y_ydy_array(const uint_fast32_t nmax,
                                           const DATA_TYPE z, DATA_TYPE *y,
                                           DATA_TYPE *dy) {

    // compute the 1st and 2nd order functions manually
    const DATA_TYPE cosz = std::cos(z);
    const DATA_TYPE sinz = std::sin(z);
    const DATA_TYPE zinv = 1. / z;
    const DATA_TYPE zinv2 = zinv * zinv;
    const DATA_TYPE zinv3 = zinv2 * zinv;
    y[0] = -cosz * zinv2 - sinz * zinv;
    y[1] = (-3. * zinv3 + zinv) * cosz - 3. * zinv2 * sinz;
    // same for the derivatives (this implicitly uses the recursion relation)
    dy[0] = -zinv * (cosz + y[0]);
    dy[1] = y[0] - zinv * y[1];
    // now apply the recursion relations for the rest
    for (uint_fast32_t i = 2; i < nmax; ++i) {
      y[i] = (2. * i + 1.) * zinv * y[i - 1] - y[i - 2];
      dy[i] = y[i - 1] - (i + 1.) * zinv * y[i];
    }
  }

  /**
   * @brief Spherical Bessel function of the first kind for general real or
   * complex numbers that returns all spherical Bessel functions and a special
   * form of their first derivatives up to the given maximum order.
   *
   * Since this method is not part of the standard C++ library, we have to
   * implement it ourselves. We make use of the following two recursion
   * relations for spherical Bessel functions:
   * @f[
   *    j_{n-1}(z) + j_{n+1}(z) = \frac{2n + 1}{z} j_n(z),
   * @f]
   * and
   * @f[
   *    \frac{d}{dz} j_n(z) = j_{n-1}(z) - \frac{n+1}{z} j_n(z),
   * @f]
   * combined with our knowledge of the zeroth order spherical Bessel function
   * of the first kind:
   * @f[
   *    j_0(z) = \frac{\sin(z)}{z}.
   * @f]
   *
   * Unfortunately, an algorithm based on a forward recursion of the zeroth
   * order function is not stable for large values of @f$n@f$ and small values
   * of @f$z@f$. Because of this, we make use of the following variant of the
   * first recursion relation:
   * @f[
   *    \frac{1}{\rho{}_n(z)} + \rho{}_{n+1}(z) = \frac{2n+1}{z},
   * @f]
   * with @f$\rho{}_n(z) = j_n(z) / j_{n-1}(z)@f$. This because the ratio of
   * two successive Bessel functions is generally more well-behaved than the
   * individual Bessel functions for the whole domain.
   *
   * We use a backward recursion algorithm for @f$\rho_n(z)@f$ that exploits
   * the following relation for large @f$n@f$ (we use @f$n = 800@f$):
   * @f[
   *    \frac{j_n(z)}{j_{n-1}(z)} = \rho_n(z) \sim{} \frac{z}{2n+1},
   * @f]
   * after which we use the recursion relation for @f$\rho{}_n(z)@f$ to find
   * the ratio for lower values of @f$n@f$.
   *
   * After this, we use a forward algorithm to determine all values of
   * @f$j_n(z)@f$ and we determine the first derivatives using the second
   * recursion relation.
   *
   * Since we do not require the actual derivative, but rather the expression
   * @f[
   *    \frac{(zj_n(z))'}{z} = j_n'(z) + \frac{j_n(z)}{z},
   * @f]
   * we make use of the fact that this expression can be computed very easily
   * during the recursion step and skip the explicit computation of the actual
   * derivative.
   *
   * @param nmax Maximum order to compute (we compute all @f$j_n(z)@f$ for
   * @f$n\in{}[1,n_{max}]@f$).
   * @param z Input value.
   * @param j Array to store the Bessel function values in (of size nmax).
   * @param dj Array to store the first derivatives in (of size nmax).
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void spherical_j_jdj_array(const uint_fast32_t nmax,
                                           const DATA_TYPE z, DATA_TYPE *j,
                                           DATA_TYPE *dj) {

    // set up and compute the array of ratios using a backward recursion
    // algorithm
    DATA_TYPE rho[BESSELFUNCTIONS_NMAX];
    const DATA_TYPE zinv = 1. / z;
    // we assume that for high enough order, the ratio tends to the
    // asymptotic value
    rho[BESSELFUNCTIONS_NMAX - 1] = z / (2. * BESSELFUNCTIONS_NMAX + 1.);
    // now recurse down to get the ratio for lower orders
    for (uint_fast32_t i = 1; i < BESSELFUNCTIONS_NMAX; ++i) {
      const uint_fast32_t index = BESSELFUNCTIONS_NMAX - i;
      rho[index - 1] = 1. / ((2. * index + 1.) * zinv - rho[index]);
    }
    // compute the zeroth order Bessel function of the first kind
    const DATA_TYPE sinz = std::sin(z);
    const DATA_TYPE j0 = sinz * zinv;
    // use the ratio and recursion relation to find the 1st order function and
    // derivative expression
    j[0] = rho[0] * j0;
    dj[0] = j0 - j[0] * zinv;
    // now recurse forward to find the rest
    for (uint_fast32_t i = 1; i < nmax; ++i) {
      j[i] = rho[i] * j[i - 1];
      dj[i] = j[i - 1] - (i + 1.) * zinv * j[i];
    }
    return;
  }
};
