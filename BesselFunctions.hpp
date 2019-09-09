#include <cinttypes>
#include <cmath>

#define BESSELFUNCTIONS_NMAX 800u

class BesselFunctions {

public:
  static inline double spherical_j(const uint_fast32_t n, const double x) {
    return std::sph_bessel(n, x);
  }

  static inline double spherical_dj(const uint_fast32_t n, const double x) {
    const double jnm1 = std::sph_bessel(n - 1, x);
    const double jn = std::sph_bessel(n, x);
    return jnm1 - (n + 1.) * jn / x;
  }

  static inline void spherical_j_dj_array(const uint_fast32_t nmax,
                                          const double x, double *j,
                                          double *dj) {
    j[0] = std::sph_bessel(1, x);
    const double xinv = 1. / x;
    dj[0] = std::sph_bessel(0, x) - 2. * xinv * j[0];
    for (uint_fast32_t i = 1; i < nmax; ++i) {
      j[i] = std::sph_bessel(i + 1, x);
      dj[i] = j[i - 1] - (i + 2.) * xinv * j[i];
    }
  }

  inline static double spherical_y(const uint_fast32_t n, const double x) {
    return std::sph_neumann(n, x);
  }

  inline static double spherical_dy(const uint_fast32_t n, const double x) {
    const double ynm1 = std::sph_neumann(n - 1, x);
    const double yn = std::sph_neumann(n, x);
    return ynm1 - (n + 1.) * yn / x;
  }

  static inline void spherical_y_dy_array(const uint_fast32_t nmax,
                                          const double x, double *y,
                                          double *dy) {
    y[0] = std::sph_neumann(1, x);
    const double xinv = 1. / x;
    dy[0] = std::sph_neumann(0, x) - 2. * xinv * y[0];
    for (uint_fast32_t i = 1; i < nmax; ++i) {
      y[i] = std::sph_neumann(i + 1, x);
      dy[i] = y[i - 1] - (i + 2.) * xinv * y[i];
    }
  }

  /**
   * @brief Spherical Bessel function of the first kind for complex numbers
   * that returns all spherical Bessel functions and their first derivatives
   * up to the given maximum order.
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
   * the following relation for large @f$n@f$ (we use @f$n = 1000@f$):
   * @f[
   *    \frac{j_n(z)}{j_{n-1}(z)} = \rho_n(z) \sim{} \frac{z}{2n+1},
   * @f]
   * after which we use the recursion relation for @$f\rho{}_n(z)@f$ to find
   * the ratio for lower values of @f$n@f$.
   *
   * After this, we use a forward algorithm to determine all values of
   * @f$j_n(z)@f$ and we determine the first derivatives using the second
   * recursion relation.
   *
   * Throughout the function, we repeatedly need to compute the inverse of
   * a complex number, @$fz = x + iy@f$, using the relation
   * @f[
   *    \frac{1}{x + iy} = \frac{x - iy}{x^2 + y^2}.
   * @f]
   * We also use
   * @f[
   *    \sin(x+iy) = \sin(x) \cosh(y) + i \cos(x) \sinh(y).
   * @f]
   *
   * @param nmax Maximum order to compute (we compute all @f$j_n(z)@f$ for
   * @f$n\in{}[1,n_{max}]@f$).
   * @param xreal Real parts of the input values (@f$x@f$ in @f$z = x + iy@f$).
   * @param ximag Imaginary parts of the input values (@f$y@f$ in
   * @f$z = x + iy@f$).
   * @param jreal Array to store the real parts of the Bessel function values
   * in (of size nmax, @f$j_r@f$ in @f$j = j_r + ij_i@f$).
   * @param jimag Array to store the imaginary parts of the Bessel function
   * values in (of size nmax, @f$j_i@f$ in @f$j = j_r + ij_i@f$).
   * @param djreal Array to store the real parts of the first derivatives in (of
   * size nmax, @f$dj_r@f$ in @f$dj = dj_r + idj_i@f$).
   * @param djreal Array to store the imaginary parts of the first derivatives
   * in (of size nmax, @f$dj_i@f$ in @f$dj = dj_r + idj_i@f$).
   */
  static inline void
  spherical_j_dj_complex_array(const uint_fast32_t nmax, const double xreal,
                               const double ximag, double *jreal, double *jimag,
                               double *djreal, double *djimag) {

    double rho_real[BESSELFUNCTIONS_NMAX];
    double rho_imag[BESSELFUNCTIONS_NMAX];

    const double zinv_denom = 1. / (xreal * xreal + ximag * ximag);
    const double zinv_real = xreal * zinv_denom;
    const double zinv_imag = -ximag * zinv_denom;
    const double proportionality_factor = 1. / (2. * BESSELFUNCTIONS_NMAX + 1.);
    rho_real[BESSELFUNCTIONS_NMAX - 1] = xreal * proportionality_factor;
    rho_imag[BESSELFUNCTIONS_NMAX - 1] = ximag * proportionality_factor;
    for (uint_fast32_t i = 1; i < BESSELFUNCTIONS_NMAX; ++i) {
      const uint_fast32_t index = BESSELFUNCTIONS_NMAX - i;
      const double recursion_factor = 2. * index + 1.;
      const double rhoinv_real = zinv_real * recursion_factor - rho_real[index];
      const double rhoinv_imag = zinv_imag * recursion_factor - rho_imag[index];
      const double rho_denom =
          1. / (rhoinv_real * rhoinv_real + rhoinv_imag * rhoinv_imag);
      rho_real[index - 1] = rhoinv_real * rho_denom;
      rho_imag[index - 1] = -rhoinv_imag * rho_denom;
    }
    //    const double sinz_real = std::sin(xreal) * std::cosh(ximag);
    //    const double sinz_imag = std::cos(xreal) * std::sinh(ximag);
    //    const double j0real = zinv_real * sinz_real - zinv_imag * sinz_imag;
    //    const double j0imag = zinv_imag * sinz_real + zinv_real * sinz_imag;
    const double rho0real_inv = zinv_real - rho_real[0];
    const double rho0imag_inv = zinv_imag - rho_imag[0];
    const double rho0_denom =
        1. / (rho0real_inv * rho0real_inv + rho0imag_inv * rho0imag_inv);
    const double rho0real = rho0real_inv * rho0_denom;
    const double rho0imag = -rho0imag_inv * rho0_denom;
    const double cosz_real = std::cos(xreal) * std::cosh(ximag);
    const double cosz_imag = -std::sin(xreal) * std::sinh(ximag);
    const double jm1real = cosz_real * zinv_real - cosz_imag * zinv_imag;
    const double jm1imag = cosz_imag * zinv_real + cosz_real * zinv_imag;
    const double j0real = rho0real * jm1real - rho0imag * jm1imag;
    const double j0imag = rho0imag * jm1real + rho0real * jm1imag;
    // remember: j[0] is actually j_1!!
    jreal[0] = rho_real[0] * j0real - rho_imag[0] * j0imag;
    jimag[0] = rho_imag[0] * j0real + rho_real[0] * j0imag;
    djreal[0] = j0real - 2. * (zinv_real * jreal[0] - zinv_imag * jimag[0]);
    djimag[0] = j0imag - 2. * (zinv_imag * jreal[0] + zinv_real * jimag[0]);
    // now recurse until nmax
    for (uint_fast32_t i = 1; i < nmax; ++i) {
      jreal[i] = rho_real[i] * jreal[i - 1] - rho_imag[i] * jimag[i - 1];
      jimag[i] = rho_imag[i] * jreal[i - 1] + rho_real[i] * jimag[i - 1];
      djreal[i] = jreal[i - 1] -
                  (i + 2.) * (zinv_real * jreal[i] - zinv_imag * jimag[i]);
      djimag[i] = jimag[i - 1] -
                  (i + 2.) * (zinv_imag * jreal[i] + zinv_real * jimag[i]);
    }
  }
};
