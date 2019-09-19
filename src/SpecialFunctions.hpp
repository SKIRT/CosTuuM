/**
 * @file SpecialFunctions.hpp
 *
 * @brief Class namespace that contains a number of special mathematical
 * functions that are used by the T-matrix code: spherical Bessel functions of
 * the first and second kind, and Wigner D functions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPECIALFUNCTIONS_HPP
#define SPECIALFUNCTIONS_HPP

#include <cinttypes>
#include <cmath>
#include <complex>
#include <vector>

/*! @brief Starting order for the backwards recurrence algorithm for the
 *  spherical Bessel functions of the first kind. The higher this number, the
 *  higher the accuracy of the Bessel functions of the first kind, but also
 *  the slower the function becomes. */
const uint_fast32_t SPECIALFUNCTIONS_BESSEL_NMAX = 800;

/**
 * @brief Class namespace that contains template functions to compute general
 * spherical Bessel functions of the first and second kind and a variant of
 * their derivatives, for real and complex input values.
 */
class SpecialFunctions {

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

    // we need these for std::complex DATA_TYPEs with extended precision,
    // since otherwise the compiler cannot find the correct division and
    // multiplication operators
    const DATA_TYPE one(1.);
    const DATA_TYPE three(3.);

    // compute the 1st and 2nd order functions manually
    const DATA_TYPE cosz = cos(z);
    const DATA_TYPE sinz = sin(z);
    const DATA_TYPE zinv = one / z;
    const DATA_TYPE zinv2 = zinv * zinv;
    const DATA_TYPE zinv3 = zinv2 * zinv;
    y[0] = -cosz * zinv2 - sinz * zinv;
    y[1] = (-three * zinv3 + zinv) * cosz - three * zinv2 * sinz;
    // same for the derivatives (this implicitly uses the recursion relation)
    dy[0] = -zinv * (cosz + y[0]);
    dy[1] = y[0] - zinv * y[1];
    // now apply the recursion relations for the rest
    for (uint_fast32_t i = 2; i < nmax; ++i) {
      const DATA_TYPE twoip1(2. * i + 1.);
      const DATA_TYPE ip1(i + 1.);
      y[i] = twoip1 * zinv * y[i - 1] - y[i - 2];
      dy[i] = y[i - 1] - ip1 * zinv * y[i];
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
    DATA_TYPE rho[SPECIALFUNCTIONS_BESSEL_NMAX];
    const DATA_TYPE one(1.);
    const DATA_TYPE zinv = one / z;
    // we assume that for high enough order, the ratio tends to the
    // asymptotic value
    const DATA_TYPE twonp1(2. * SPECIALFUNCTIONS_BESSEL_NMAX + 1.);
    rho[SPECIALFUNCTIONS_BESSEL_NMAX - 1] = z / twonp1;
    // now recurse down to get the ratio for lower orders
    for (uint_fast32_t i = 1; i < SPECIALFUNCTIONS_BESSEL_NMAX; ++i) {
      const uint_fast32_t index = SPECIALFUNCTIONS_BESSEL_NMAX - i;
      const DATA_TYPE twoip1(2. * index + 1.);
      rho[index - 1] = one / (twoip1 * zinv - rho[index]);
    }
    // compute the zeroth order Bessel function of the first kind
    const DATA_TYPE sinz = sin(z);
    const DATA_TYPE j0 = sinz * zinv;
    // use the ratio and recursion relation to find the 1st order function and
    // derivative expression
    j[0] = rho[0] * j0;
    dj[0] = j0 - j[0] * zinv;
    // now recurse forward to find the rest
    for (uint_fast32_t i = 1; i < nmax; ++i) {
      j[i] = rho[i] * j[i - 1];
      const DATA_TYPE ip1(i + 1.);
      dj[i] = j[i - 1] - ip1 * zinv * j[i];
    }
    return;
  }

  /**
   * @brief Compute the Wigner D functions and derivatives for the given quantum
   * number @f$m@f$ and up to the given order @f$n_{max}@f$.
   *
   * This function returns
   * @f[
   *    d^n_{0m}(x) = (-1)^{-m} \sqrt{\frac{(n-m)!}{(n+m)!}} P^m_n(\cos(x)),
   * @f]
   * where @f$n \in{} [1, n_{max}]@f$, @f$m \in{} [-n, n]@f$ and
   * @f$P^m_n(x)@f$ is the associated Legendre polynomial of degree @f$n@f$ and
   * order @f$m@f$ (see
   * https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
   *
   * We compute the function using a forward recursion algorithm:
   * @f[
   *    d^{n+1}_{0m}(x) = \frac{1}{\sqrt{(n+1)^2 - m^2}} \left(
   *      (2n + 1) \cos(x) d^n_{0m}(x) - \sqrt{n^2 - m^2} d^{n-1}_{0m}(x)
   *      \right)
   * @f]
   * and
   * @f[
   *    \frac{d}{dx} d^n_{0m}(x) = \frac{1}{(2n+1)\sin(x)} \left(
   *      -(n+1)\sqrt{n^2 - m^2} d^{n-1}_{0m}(x)
   *      + n \sqrt{(n+1)^2 - m^2} d^{n+1}_{0m}(x) \right).
   * @f]
   * As initial conditions, we use
   * @f[
   *    d^{m-1}_{0m} (x) = 0
   * @f]
   * and
   * @f[
   *    d^m_{0m} (x) = A_m \sin^m(x),
   * @f]
   * with
   * @f[
   *    A_{m+1} = A_m \sqrt{\frac{2m + 1}{2(m+1)}},
   * @f]
   * and @f$A_0 = 1@f$. These equations can be found in Mishchenko, M. I., 2000,
   * Applied Optics, 39, 1026 (https://doi.org/10.1364/AO.39.001026). They can
   * also be derived from the recursion relations for the associated Legendre
   * polynomials that can be found on Wikipedia. Note that there is a missing
   * @f$n@f$ factor in the second term of the second recursion relation in the
   * Mishchenko paper.
   *
   * Note that we set all @f$d^{m}_{0n}(x)@f$ with @f$n < m-1@f$ equal to zero
   * in the return arrays. Also note that the code for @f$m = 0@f$ differs from
   * the code for general @f$m > 0@f$, since we require a different starting
   * point for the recursion algorithm.
   *
   * @param cosx Cosine of the input value @f$x@f$ (should be in the interval
   * @f$[-1,1]@f$).
   * @param nmax Maximum order @f$n_{max}@f$.
   * @param m Absolute value of the quantum number @f$m@f$
   * (@f$m \in{} [0,n_{max}]@f$).
   * @param y Array to store the functions @f$d^{n}_{0m}(x)@f$ in (of size
   * nmax).
   * @param dy Array to store the derivatives
   * @f$d/dx\left(d^n_{0m}(x)\right)@f$ in (of size nmax).
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void
  wigner_dn_0m(const DATA_TYPE cosx, const uint_fast32_t nmax,
               const uint_fast32_t m, DATA_TYPE *y, DATA_TYPE *dy) {

    const DATA_TYPE zero(0.);
    const DATA_TYPE one(1.);
    const DATA_TYPE two(2.);

    // precompute sin(x) and its inverse
    const DATA_TYPE sinx = sqrt(one - cosx * cosx);
    const DATA_TYPE sinxinv = one / sinx;
    // branch out depending on the m value
    if (m == 0) {
      // manually set the starting point for the recursion algorithm by
      // applying the recursion formula once
      DATA_TYPE d1(1.);
      DATA_TYPE d2 = cosx;
      // now recurse for the other values
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        const DATA_TYPE np1(n + 1.);
        const DATA_TYPE np2(n + 2.);
        const DATA_TYPE np12p1 = two * np1 + one;
        const DATA_TYPE dnp1 = (np12p1 * cosx * d2 - np1 * d1) / np2;
        const DATA_TYPE ddn = sinxinv * np1 * np2 * (dnp1 - d1) / np12p1;
        // note that we computed d^{n+1}_{0m}; d^n_{0m} was computed in the
        // previous iteration
        y[n] = d2;
        dy[n] = ddn;
        d1 = d2;
        d2 = dnp1;
      }
    } else {
      // precompute m^2
      const DATA_TYPE m2 = m * m;
      DATA_TYPE d1(0.);
      DATA_TYPE d2(1.);
      // compute the factor A_m sin^m(x)
      for (uint_fast32_t i = 0; i < m; ++i) {
        const DATA_TYPE i2(2. * (i + 1.));
        d2 *= sqrt((i2 - one) / i2) * sinx;
      }
      // zero out n < m-1 elements in the array
      for (uint_fast32_t n = 0; n < m - 1; ++n) {
        y[n] = zero;
        dy[n] = zero;
      }
      // now recurse to find all other values
      for (uint_fast32_t n = m - 1; n < nmax; ++n) {
        const DATA_TYPE np1(n + 1.);
        const DATA_TYPE np2(n + 2.);
        const DATA_TYPE np12p1 = two * np1 + one;
        const DATA_TYPE sqrtnp12mm2 = sqrt(np1 * np1 - m2);
        const DATA_TYPE sqrtnp22mm2 = sqrt(np2 * np2 - m2);
        const DATA_TYPE dnp1 =
            (np12p1 * cosx * d2 - sqrtnp12mm2 * d1) / sqrtnp22mm2;
        const DATA_TYPE ddn =
            sinxinv * (np1 * sqrtnp22mm2 * dnp1 - np2 * sqrtnp12mm2 * d1) /
            np12p1;
        // note that we computed d^{n+1}_{0m}; d^n_{0m} was computed in the
        // previous iteration
        y[n] = d2;
        dy[n] = ddn;
        d1 = d2;
        d2 = dnp1;
      }
    }
  }

  /**
   * @brief Special version of SpecialFunctions::wigner_dn_0m() that divides the
   * elements of the first return array by the sine of the input angle.
   *
   * The code is mostly the same as the original function with the addition of
   * the division by the sine. However, for input angles close to zero, possible
   * division by zero and serious round off error might occur, which we need to
   * deal with here with special code.
   *
   * @param cosx Cosine of the input value @f$x@f$ (should be in the interval
   * @f$[-1,1]@f$).
   * @param nmax Maximum order @f$n_{max}@f$.
   * @param m Absolute value of the quantum number @f$m@f$
   * (@f$m \in{} [0,n_{max}]@f$).
   * @param y Array to store the functions @f$\frac{d^{n}_{0m}(x)}{\sin(x)}@f$
   * in (of size nmax).
   * @param dy Array to store the derivatives
   * @f$d/dx\left(d^n_{0m}(x)\right)@f$ in (of size nmax).
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void
  wigner_dn_0m_sinx(const DATA_TYPE cosx, const uint_fast32_t nmax,
                    const uint_fast32_t m, DATA_TYPE *y, DATA_TYPE *dy) {

    // check if the input argument is close to 1 (corresponding to a sine close
    // to 0
    if (abs(1. - abs(cosx)) <= 1.e-10) {
      if (m == 1) {
        int_fast8_t sign = 1;
        for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
          DATA_TYPE dn = 0.5 * sqrt(n * (n + 1.));
          if (cosx < 0.) {
            dn *= sign;
          }
          y[n - 1] = dn;
          if (cosx < 0.) {
            dn = -dn;
          }
          dy[n - 1] = dn;
          // swap sign for the next iteration
          sign = -sign;
        }
      } else {
        // simply set all return values to zero
        for (uint_fast32_t i = 0; i < nmax; ++i) {
          y[i] = 0.;
          dy[i] = 0.;
        }
      }
    } else {
      const DATA_TYPE zero(0.);
      const DATA_TYPE one(1.);
      const DATA_TYPE two(2.);

      // precompute sin(x) and its inverse
      const DATA_TYPE sinx = sqrt(one - cosx * cosx);
      const DATA_TYPE sinxinv = one / sinx;
      // branch out depending on the m value
      if (m == 0) {
        // manually set the starting point for the recursion algorithm by
        // applying the recursion formula once
        DATA_TYPE d1(1.);
        DATA_TYPE d2 = cosx;
        // now recurse for the other values
        for (uint_fast32_t n = 0; n < nmax; ++n) {
          const DATA_TYPE np1(n + 1.);
          const DATA_TYPE np2(n + 2.);
          const DATA_TYPE np12p1 = two * np1 + one;
          const DATA_TYPE dnp1 = (np12p1 * cosx * d2 - np1 * d1) / np2;
          const DATA_TYPE ddn = sinxinv * np1 * np2 * (dnp1 - d1) / np12p1;
          // note that we computed d^{n+1}_{0m}; d^n_{0m} was computed in the
          // previous iteration
          y[n] = d2 * sinxinv;
          dy[n] = ddn;
          d1 = d2;
          d2 = dnp1;
        }
      } else {
        // precompute m^2
        const DATA_TYPE m2 = m * m;
        DATA_TYPE d1(0.);
        DATA_TYPE d2(1.);
        // compute the factor A_m sin^m(x)
        for (uint_fast32_t i = 0; i < m; ++i) {
          const DATA_TYPE i2(2. * (i + 1.));
          d2 *= sqrt((i2 - one) / i2) * sinx;
        }
        // zero out n < m-1 elements in the array
        for (uint_fast32_t n = 0; n < m - 1; ++n) {
          y[n] = zero;
          dy[n] = zero;
        }
        // now recurse to find all other values
        for (uint_fast32_t n = m - 1; n < nmax; ++n) {
          const DATA_TYPE np1(n + 1.);
          const DATA_TYPE np2(n + 2.);
          const DATA_TYPE np12p1 = two * np1 + one;
          const DATA_TYPE sqrtnp12mm2 = sqrt(np1 * np1 - m2);
          const DATA_TYPE sqrtnp22mm2 = sqrt(np2 * np2 - m2);
          const DATA_TYPE dnp1 =
              (np12p1 * cosx * d2 - sqrtnp12mm2 * d1) / sqrtnp22mm2;
          const DATA_TYPE ddn =
              sinxinv * (np1 * sqrtnp22mm2 * dnp1 - np2 * sqrtnp12mm2 * d1) /
              np12p1;
          // note that we computed d^{n+1}_{0m}; d^n_{0m} was computed in the
          // previous iteration
          y[n] = d2 * sinxinv;
          dy[n] = ddn;
          d1 = d2;
          d2 = dnp1;
        }
      }
    }
  }

  /**
   * @brief Calculates the ratio of the radii of an equal volume and an equal
   * surface area sphere for the spheroid with the given axis ratio.
   *
   * Assume a general spheroid with equation
   * @f[
   *    \frac{x^2}{a^2} + \frac{y^2}{a^2} + \frac{z^2}{b^2} = 1.
   * @f]
   * If @f$a > b@f$, we call this spheroid oblate. In this case, the
   * eccentricity is defined as @f$e^2 = 1 - \frac{b^2}{a^2}@f$, and the surface
   * area of the spheroid is given by
   * @f[
   *    S_{obl} = 2\pi{} a^2
   *              + \frac{\pi{}b^2}{e}\ln\left(\frac{1+e}{1-e}\right).
   * @f]
   * If @f$a < b@f$, the spheroid is prolate, its eccentricity is
   * @f$e^2 = 1 - \frac{a^2}{b^2}@f$, and the surface area is given by
   * @f[
   *    S_{prol} = 2\pi{}a^2 + \frac{2\pi{}ab}{e}\arcsin(e).
   * @f]
   * In both cases, the volume of the spheroid is given by
   * @f[
   *    V_{sph} = \frac{4\pi{}}{3}a^2 b.
   * @f]
   *
   * We want to find the two equivalent sphere radii @f$R_V@f$ and @f$R_S@f$,
   * defined by imposing that either @f$V = \frac{4\pi{}}{3}R_V^3@f$ or
   * @f$S = 4\pi{}R_S^2@f$ equals its spheroid counterpart. This means that the
   * equivalent sphere radius @f$R_V@f$ is simply given by
   * @f[
   *    R_V = a^{\frac{2}{3}}b^{\frac{1}{3}}.
   * @f]
   * Depending on the type of spheroid, the equivalent sphere radius @f$R_S@f$
   * is either
   * @f[
   *    R_{S,obl} = \sqrt{\frac{a^2}{2}
   *                + \frac{b^2}{4e}\ln\left(\frac{1+e}{1-e}\right)},
   * @f]
   * or
   * @f[
   *    R_{S,prol} = \sqrt{\frac{a^2}{2} + \frac{ab}{2e}\arcsin(e)}.
   * @f]
   * The ratio @f$\frac{R_V}{R_S}@f$ is then either
   * @f[
   *    \frac{R_V}{R_{S,obl}} = \left(\sqrt{\frac{d^{\frac{2}{3}}}{2}
   *  + \frac{d^{-\frac{4}{3}}}{4e}\ln\left(\frac{1+e}{1-e}\right)}\right)^{-1},
   * @f]
   * or
   * @f[
   *    \frac{R_V}{R_{S,prol}} = \left(\sqrt{\frac{d^{\frac{2}{3}}}{2}
   *  + \frac{d^{-\frac{1}{3}}}{2e}\arcsin(e)}\right)^{-1},
   * @f]
   * where we substituted @f$d = \frac{a}{b}@f$.
   *
   * @param axis_ratio Ratio of the horizontal and vertical axis of the
   * spheroid, @f$d = \frac{a}{b}@f$.
   * @return Ratio of the radii of the equal volume and equal surface area
   * sphere, @f$\frac{R_V}{R_S}@f$.
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline DATA_TYPE get_equal_volume_to_equal_surface_area_sphere_ratio(
      const DATA_TYPE axis_ratio) {

    const DATA_TYPE one(1.);
    const DATA_TYPE two(2.);
    const DATA_TYPE half(0.5);
    const DATA_TYPE quarter(0.25);

    if (axis_ratio > one) {
      const DATA_TYPE axis_ratio2 = axis_ratio * axis_ratio;
      const DATA_TYPE e = sqrt(one - one / axis_ratio2);
      const DATA_TYPE cbrtaxis_ratio2 = cbrt(axis_ratio2);
      return one /
             sqrt(quarter * (two * cbrtaxis_ratio2 +
                             log((one + e) / (one - e)) /
                                 (e * cbrtaxis_ratio2 * cbrtaxis_ratio2)));
    } else if (axis_ratio < one) {
      const DATA_TYPE cbrtaxis_ratio = cbrt(axis_ratio);
      const DATA_TYPE e = sqrt(one - axis_ratio * axis_ratio);
      return one / sqrt(half * (cbrtaxis_ratio * cbrtaxis_ratio +
                                asin(e) / (e * cbrtaxis_ratio)));
    } else {
      // the spheroid is in fact as sphere, so the ratio is trivially 1
      return one;
    }
  }

  /**
   * @brief Get the coordinates and weights for a 1D Gauss-Legendre quadrature
   * on the interval @f$[-1,1]@f$ with the given order.
   *
   * We are looking for the roots @f$x_i@f$ of the Legendre polynomial of order
   * @f$n@f$, defined via the recursion relation
   * @f[
   *    P_n(x) = \frac{2n-1}{n}xP_{n-1}(x) - \frac{n-1}{n}P_{n-2}(x)
   *           = xP_{n-1}(x) + \frac{n-1}{n}
   *                \left(xP_{n-1}(x) - P_{n-2}(x)\right),
   * @f]
   * with
   * @f[
   *    P_0(x) = 1
   * @f]
   * and
   * @f[
   *    P_1(x) = x.
   * @f]
   * For each of these roots, the associated Gauss-Legendre weight is given by
   * @f[
   *    w_i = \frac{2}{(1-x_i^2)P_n'^2(x_i)},
   * @f]
   * where the derivative of the Legendre polynomial of order @f$n@f$ satisfies
   * @f[
   *    P_n'(x) = \frac{n}{x^2-1}\left(xP_n(x) - P_{n-1}(x)\right).
   * @f]
   *
   * The weights are chosen such that an approximate quadrature rule
   * @f[
   *    \int_{-1}^{1} f(x) dx \approx{} \sum_{i=1}^{n} w_i f(x_i),
   * @f]
   * is exact for polynomials of up to degree @f$2n-1@f$, this is called
   * Gauss-Legendre quadrature.
   *
   * To compute the roots, we use a Newton-Raphson iterative scheme. We start
   * with a suitable initial guess for each zero point @f$x_i@f$, and then use
   * the two recursion relations to update this initial guess:
   * @f[
   *    x_i \rightarrow{} x_i - \frac{P_n(x_i)}{P_n'{x_i}}
   * @f]
   * until the value @f$P_n(x_i)@f$ becomes very small. As initial guess, we
   * use
   * @f[
   *    x_i = \left[1 - \frac{1}{8}\frac{1}{n^2} + \frac{5}{38}\frac{1}{n^3}
   *                - \frac{2}{25}\frac{1}{n^4}
   *                    \left(1 - \frac{14}{39}\frac{1}{\theta_{n,i}^2} \right)
   *          \right] \cos\left( \theta_{n,i} \right),
   * @f]
   * where @f$i = 1...n@f$ and @f$\theta_{n,i}=\pi{}
   * \frac{i-\frac{1}{4}}{n + \frac{1}{2}}@f$; this expression is based on
   * Lether & Wenston, 1995, Journal of Computational and Applied Mathematics,
   * 59, 245 (https://doi.org/10.1016/0377-0427(94)00030-5).
   *
   * @param order Order of the quadrature, @f$n@f$.
   * @param points Vector to store the resulting coordinates in. This vector
   * should be preallocated with the correct size, i.e. the order @f$n@f$.
   * @param weights Vector to store the resulting weights in. This vector
   * should be preallocated with the correct size, i.e. the order @f$n@f$.
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void
  get_gauss_legendre_points_and_weights(const uint_fast32_t order,
                                        std::vector<DATA_TYPE> &points,
                                        std::vector<DATA_TYPE> &weights) {

    // we know that the roots are symmetric around 0 in the interval [-1,1]
    // this means we only need to compute half of them (+1 for odd n)
    const uint_fast32_t ind = order % 2;
    const uint_fast32_t k = order / 2 + ind;
    for (uint_fast32_t i = 0; i < k; ++i) {

      // compute the initial guess for x_i (note that our i is i-1)
      const DATA_TYPE theta_n_k = M_PI * (i + 0.75) / (order + 0.5);
      const DATA_TYPE orderinv = 1. / order;
      const DATA_TYPE orderinv2 = orderinv * orderinv;
      DATA_TYPE x_i =
          (1. - 0.125 * orderinv2 + (5. / 38.) * orderinv * orderinv2 -
           (2. / 25.) * orderinv2 * orderinv2 *
               (1. - (14. / 39.) / (theta_n_k * theta_n_k))) *
          cos(theta_n_k);

      // Newton-Raphson scheme
      uint_fast32_t niter = 0;
      // we aim for a precision of 1.e-16, and relax this if we cannot obtain
      // an accurate result after many iterations
      DATA_TYPE check = 1.e-16;
      DATA_TYPE pa;
      DATA_TYPE pb;
      do {
        ++niter;
        // relax the tolerance after many iterations
        if (niter > 100) {
          check *= 10.;
        }
        // initial values for the recursion relation
        pb = 1.;
        DATA_TYPE pc = x_i;
        // now recurse until we reach the desired order
        for (uint_fast32_t j = 1; j < order; ++j) {
          pa = pb;
          pb = pc;
          pc = x_i * pb + (x_i * pb - pa) * j / (j + 1.);
        }
        // pc now contains P_n(x)
        pa = 1. / ((pb - x_i * pc) * order);
        // pa now contains 1/[(1-x^2)P_n'(x)]
        pb = pa * pc * (1. - x_i * x_i);
        // pb is P_n(x)/P_n'(x)
        x_i = x_i - pb;
      } while (std::abs(pb) > check * std::abs(x_i));

      // we set the element n - i - 1
      const uint_fast32_t m = order - i - 1;
      points[m] = x_i;
      weights[m] = 2. * pa * pa * (1. - x_i * x_i);
      // if this is not the final root in an odd n calculation, we also set
      // element i
      if (i != k - 1 || ind != 1) {
        points[i] = -points[m];
        weights[i] = weights[m];
      }
    }
  }

  /**
   * @brief Get the radius (squared) and derivative w.r.t. azimuthal angle
   * divided by the radius for an spheroid with the given equal volume sphere
   * radius and axis ratio, for the given input azimuthal angles.
   *
   * The equation of a spheroid is
   * @f[
   *    \frac{x^2}{a^2} + \frac{y^2}{a^2} + \frac{z^2}{b^2} = 1.
   * @f]
   * If we substitute the axis ratio @f$d = \frac{a}{b}@f$ and use the
   * expressions for @f$x, y, z@f$ in spherical coordinates, we get the
   * following equation for the spheroid:
   * @f[
   *    r^2 = \frac{a^2}{\sin^2(\theta{}) + d^2 \cos^2(\theta{})}.
   * @f]
   * The horizontal axis @f$a@f$ can be found from the definition of the equal
   * volume sphere radius:
   * @f[
   *    V_S = \frac{4\pi{}}{3}R_V^2 = \frac{4\pi{}}{3}a^2b
   *        = \frac{4\pi{}}{3} \frac{a^3}{d} = V_{sph},
   * @f]
   * hence
   * @f[
   *    a = R_V d^{\frac{1}{3}}.
   * @f]
   *
   * If we take the square root of the expression for @f$r^2@f$ and derive
   * w.r.t. the azimuthal angle @f$\theta{}@f$, we get
   * @f[
   *    \frac{dr}{d\theta{}} = \frac{d}{d\theta{}} \left(
   *        \frac{a}{\sqrt{\sin^2(\theta{}) + d^2\cos^2(\theta{})}}\right)
   *      = \frac{a\sin(\theta{})\cos(\theta{})(d^2-1)}{
   *        \left(\sin^2(\theta{}) + d^2\cos^2(\theta{})\right)^{\frac{3}{2}}}
   *      = r \frac{\sin(\theta{})\cos(\theta{})(d^2-1)}{
   *        \sin^2(\theta{}) + d^2\cos^2(\theta{})}
   * @f]
   *
   * @param costheta Cosine of the azimuthal angles, @f$\cos(\theta{})@f$. We
   * assume that values are given in order from low to high, and are symmetric
   * w.r.t. @f$0@f$.
   * @param R_V Equal volume sphere radius, @f$R_V@f$.
   * @param axis_ratio Ratio of horizontal to vertical axis,
   * @f$d = \frac{a}{b}@f$.
   * @param r2 Output radii squared, @f$r^2(\theta{})@f$. This vector should be
   * preallocated with the correct size, i.e. the same size as the input vector.
   * @param dr_over_r Derivative over radius,
   * @f$\frac{1}{r(\theta{})}\frac{dr(\theta{})}{d\theta{}}@f$. This vector
   * should be preallocated with the correct size, i.e. the same size as the
   * input vector.
   * @tparam DATA_TYPE Data type of input and output values.
   */
  template <typename DATA_TYPE>
  static inline void
  get_r_dr_spheroid(const std::vector<DATA_TYPE> &costheta, const DATA_TYPE R_V,
                    const DATA_TYPE axis_ratio, std::vector<DATA_TYPE> &r2,
                    std::vector<DATA_TYPE> &dr_over_r) {

    // compute the horizontal axis length
    const DATA_TYPE a = R_V * cbrt(axis_ratio);
    const DATA_TYPE a2 = a * a;
    const DATA_TYPE axis_ratio2 = axis_ratio * axis_ratio;
    const DATA_TYPE axis_ratio2m1 = axis_ratio2 - 1.;
    for (uint_fast32_t i = 0; i < costheta.size() / 2; ++i) {
      const DATA_TYPE costheta2 = costheta[i] * costheta[i];
      const DATA_TYPE sintheta2 = 1. - costheta2;
      const DATA_TYPE sintheta = sqrt(sintheta2);
      const DATA_TYPE r2_over_a2 = 1. / (sintheta2 + axis_ratio2 * costheta2);
      r2[i] = a2 * r2_over_a2;
      r2[costheta.size() - i - 1] = r2[i];
      dr_over_r[i] = r2_over_a2 * costheta[i] * sintheta * axis_ratio2m1;
      dr_over_r[costheta.size() - i - 1] = -dr_over_r[i];
    }
  }
};

#endif // SPECIALFUNCTIONS_HPP
