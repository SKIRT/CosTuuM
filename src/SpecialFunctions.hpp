/**
 * @file SpecialFunctions.hpp
 *
 * @brief Class namespace that contains a number of special mathematical
 * functions that are used by the T-matrix code: spherical Bessel functions of
 * the first and second kind, and Wigner D functions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include <cinttypes>
#include <cmath>
#include <complex>
#include <vector>

/*! @brief Starting order for the backwards recurrence algorithm for the
 *  spherical Bessel functions of the first kind. The higher this number, the
 *  higher the accuracy of the Bessel functions of the first kind, but also
 *  the slower the function becomes. */
#define SPECIALFUNCTIONS_BESSEL_NMAX 800u

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
    DATA_TYPE rho[SPECIALFUNCTIONS_BESSEL_NMAX];
    const DATA_TYPE zinv = 1. / z;
    // we assume that for high enough order, the ratio tends to the
    // asymptotic value
    rho[SPECIALFUNCTIONS_BESSEL_NMAX - 1] =
        z / (2. * SPECIALFUNCTIONS_BESSEL_NMAX + 1.);
    // now recurse down to get the ratio for lower orders
    for (uint_fast32_t i = 1; i < SPECIALFUNCTIONS_BESSEL_NMAX; ++i) {
      const uint_fast32_t index = SPECIALFUNCTIONS_BESSEL_NMAX - i;
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
   */
  static inline void wigner_dn_0m(const double cosx, const uint_fast32_t nmax,
                                  const uint_fast32_t m, double *y,
                                  double *dy) {

    // precompute sin(x) and its inverse
    const double sinx = std::sqrt(1. - cosx * cosx);
    const double sinxinv = 1. / sinx;
    // branch out depending on the m value
    if (m == 0) {
      // manually set the starting point for the recursion algorithm by
      // applying the recursion formula once
      double d1 = 1.;
      double d2 = cosx;
      // now recurse for the other values
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        const double np1 = n + 1.;
        const double np2 = n + 2.;
        const double np12p1 = 2. * np1 + 1.;
        const double dnp1 = (np12p1 * cosx * d2 - np1 * d1) / np2;
        const double ddn = sinxinv * np1 * np2 * (dnp1 - d1) / np12p1;
        // note that we computed d^{n+1}_{0m}; d^n_{0m} was computed in the
        // previous iteration
        y[n] = d2;
        dy[n] = ddn;
        d1 = d2;
        d2 = dnp1;
      }
    } else {
      // precompute m^2
      const double m2 = m * m;
      double d1 = 0.;
      double d2 = 1.;
      // compute the factor A_m sin^m(x)
      for (uint_fast32_t i = 0; i < m; ++i) {
        const double i2 = 2. * (i + 1.);
        d2 *= std::sqrt((i2 - 1.) / i2) * sinx;
      }
      // zero out n < m-1 elements in the array
      for (uint_fast32_t n = 0; n < m - 1; ++n) {
        y[n] = 0.;
        dy[n] = 0.;
      }
      // now recurse to find all other values
      for (uint_fast32_t n = m - 1; n < nmax; ++n) {
        const double np1 = n + 1.;
        const double np2 = n + 2.;
        const double np12p1 = 2. * np1 + 1.;
        const double sqrtnp12mm2 = std::sqrt(np1 * np1 - m2);
        const double sqrtnp22mm2 = std::sqrt(np2 * np2 - m2);
        const double dnp1 =
            (np12p1 * cosx * d2 - sqrtnp12mm2 * d1) / sqrtnp22mm2;
        const double ddn = sinxinv *
                           (np1 * sqrtnp22mm2 * dnp1 - np2 * sqrtnp12mm2 * d1) /
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
   */
  static inline double
  get_equal_volume_to_equal_surface_area_sphere_ratio(const double axis_ratio) {
    if (axis_ratio > 1.) {
      const double axis_ratio2 = axis_ratio * axis_ratio;
      const double e = std::sqrt(1. - 1. / axis_ratio2);
      const double cbrtaxis_ratio2 = std::cbrt(axis_ratio2);
      return 1. /
             std::sqrt(0.25 * (2. * cbrtaxis_ratio2 +
                               std::log((1. + e) / (1. - e)) /
                                   (e * cbrtaxis_ratio2 * cbrtaxis_ratio2)));
    } else if (axis_ratio < 1.) {
      const double cbrtaxis_ratio = std::cbrt(axis_ratio);
      const double e = std::sqrt(1. - axis_ratio * axis_ratio);
      return 1. / std::sqrt(0.5 * (cbrtaxis_ratio * cbrtaxis_ratio +
                                   std::asin(e) / (e * cbrtaxis_ratio)));
    } else {
      // the spheroid is in fact as sphere, so the ratio is trivially 1
      return 1.;
    }
  }

  /**
   * @brief Get the coordinates and weights for a 1D Gauss-Legendre quadrature
   * on the interval @f$[-1,1]@f$.
   *
   * Needs to be properly documented!
   *
   * @param order Order of the quadrature.
   * @param points Vector to store the resulting coordinates in (of size order).
   * @param weights Vector to store the resulting weights in (of size order).
   */
  static inline void
  get_gauss_legendre_points_and_weigths(const uint_fast32_t order,
                                        std::vector<double> &points,
                                        std::vector<double> &weights) {
    const uint_fast32_t ind = order % 2;
    const uint_fast32_t k = order / 2 + ind;
    for (uint_fast32_t i = 0; i < k; ++i) {
      const uint_fast32_t m = order - i - 1;
      double x;
      if (i == 0) {
        x = 1. - 2. / ((order + 1.) * order);
      } else if (i == 1) {
        x = 4. * (points[order - 1] - 1.) + points[order - 1];
      } else if (i == 2) {
        x = 1.6 * (points[order - 2] - points[order - 1]) + points[order - 2];
      } else if (i == k - 1 && ind == 1) {
        x = 0.;
      } else {
        x = 3. * (points[m + 1] - points[m + 2]) + points[m + 3];
      }
      uint_fast32_t niter = 0;
      double check = 1.e-16;
      double pa;
      double pb;
      do {
        pb = 1.;
        ++niter;
        if (niter > 100) {
          check *= 10.;
        }
        double pc = x;
        for (uint_fast32_t j = 1; j < order; ++j) {
          pa = pb;
          pb = pc;
          pc = x * pb + (x * pb - pa) * j / (j + 1.);
        }
        pa = 1. / ((pb - x * pc) * order);
        pb = pa * pc * (1. - x * x);
        x = x - pb;
      } while (std::abs(pb) > check * std::abs(x));
      points[m] = x;
      weights[m] = 2. * pa * pa * (1. - x * x);
      if (i != k - 1 || ind != 1) {
        points[i] = -points[m];
        weights[i] = weights[m];
      }
    }
  }
};
