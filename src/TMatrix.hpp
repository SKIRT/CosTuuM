/**
 * @file TMatrix.hpp
 *
 * @brief T-matrix class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TMATRIX_HPP
#define TMATRIX_HPP

#include "Matrix.hpp"
#include "SpecialFunctions.hpp"

#include <complex>
#include <vector>

/**
 * @brief T-matrix and auxiliary variables required to compute it.
 */
class TMatrix {
private:
  /*! @brief Maximum order of the spherical basis functions, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of components in one quarter of the T-matrix,
   *  @f$L_{max} = n_{max} (n_{max} + 2) @f$. */
  const uint_fast32_t _Lmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Precomputed factors @f$n(n+1)@f$ (array of size @f$n_{max}@f$). */
  std::vector<double> _an;

  /*! @brief Precomputed factors @f$\sqrt{\frac{2n+1}{n(n+1)}}@f$ (array of size
   *  @f$n_{max}@f$). */
  std::vector<double> _dd;

  /*! @brief Precomputed factors
   *  @f$\frac{1}{2}\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)n'(n'+1)}}@f$
   *  (@f$n_{max}\times{}n_{max}@f$ matrix). */
  Matrix<double> _ann;

  /*! @brief Precomputed factors @f$\cos(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<double> _costheta;

  /*! @brief Precomputed factors @f$\frac{1}{\sin(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<double> _sinthetainv;

  /*! @brief Precomputed factors @f$\frac{1}{\sin^2(\theta{})}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<double> _sintheta2inv;

  /*! @brief Gauss-Legendre weights for the roots @f$\cos(\theta{})@f$ (array of
   *  size @f$2n_{GL}@f$). */
  std::vector<double> _weights;

  /*! @brief Precomputed factors @f$r^2(\theta{})@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<double> _r2;

  /*! @brief Precomputed factors @f$\frac{1}{r(\theta{})}\frac{d}{d\theta{}}
   *  r(\theta{})@f$  (array of size @f$2n_{GL}@f$). */
  std::vector<double> _dr_over_r;

  /*! @brief Precomputed factors @f$kr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<double> _kr;

  /*! @brief Precomputed factors @f$\frac{1}{kr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<double> _krinv;

  /*! @brief Precomputed factors @f$km_rr@f$ (array of size @f$2n_{GL}@f$). */
  std::vector<std::complex<double>> _krmr;

  /*! @brief Precomputed factors @f$\frac{1}{km_rr}@f$ (array of size
   *  @f$2n_{GL}@f$). */
  std::vector<std::complex<double>> _krmrinv;

  /*! @brief Wavenumber, @f$k = \frac{2\pi{}}{\lambda{}}@f$. */
  const double _k;

  /*! @brief Wavenumber squared, @f$k^2@f$. */
  const double _k2;

  /*! @brief Wavenumber squared times refractive index, @f$m_rk^2@f$. */
  const std::complex<double> _k2mr;

  /*! @brief Bessel functions @f$j_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<double> _jkr;

  /*! @brief Bessel functions @f$y_n(kr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<double> _ykr;

  /*! @brief Bessel function derivatives @f$\frac{[krj_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<double> _djkr;

  /*! @brief Bessel function derivatives @f$\frac{[kry_n(kr)]'}{kr}@f$
   *  (@f$2n_{GL}\times{}n_{max}@f$ matrix). */
  Matrix<double> _dykr;

  /*! @brief Bessel functions @f$j_n(km_rr)@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<double>> _jkrmr;

  /*! @brief Bessel function derivatives
   *  @f$\frac{[km_rrj(km_rr)]'}{km_rr}@f$ (@f$2n_{GL}\times{}n_{max}@f$
   *  matrix). */
  Matrix<std::complex<double>> _djkrmr;

  /*! @brief T-matrix itself
   *  (@f$2n_{max}(n_{max}+2)\times{}2n_{max}(n_{max}+2)@f$ matrix). */
  Matrix<std::complex<double>> _T;

public:
  /**
   * @brief T-matrix.
   *
   * @param wavelength Wavelength of incident radiation, @f$\lambda{}@f$.
   * @param refractive_index Refractive index of the material, @f$m_r@f$.
   * @param R_V Equal volume sphere radius, @f$R_V@f$.
   * @param axis_ratio Axis ratio of the spheroid, @f$d = \frac{a}{b}@f$.
   * @param nmax Maximum order of spherical basis, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline TMatrix(const double wavelength,
                 const std::complex<double> refractive_index, const double R_V,
                 const double axis_ratio, const uint_fast32_t nmax,
                 const uint_fast32_t ngauss)
      : _nmax(nmax), _Lmax(nmax * (nmax + 2)), _ngauss(ngauss), _an(nmax, 0.),
        _dd(nmax, 0.), _ann(nmax, nmax), _costheta(2 * ngauss, 0.),
        _sinthetainv(2 * ngauss, 0.), _sintheta2inv(2 * ngauss, 0.),
        _weights(2 * ngauss, 0.), _r2(2 * ngauss, 0.),
        _dr_over_r(2 * ngauss, 0.), _kr(2 * ngauss, 0.), _krinv(2 * ngauss, 0.),
        _krmr(2 * ngauss, 0.), _krmrinv(2 * ngauss, 0.),
        _k(2. * M_PI / wavelength), _k2(_k * _k), _k2mr(refractive_index * _k2),
        _jkr(2 * ngauss, nmax), _ykr(2 * ngauss, nmax), _djkr(2 * ngauss, nmax),
        _dykr(2 * ngauss, nmax), _jkrmr(2 * ngauss, nmax),
        _djkrmr(2 * ngauss, nmax), _T(2 * _Lmax, 2 * _Lmax) {

    for (uint_fast32_t ni = 0; ni < nmax; ++ni) {
      const double nn = (ni + 2.) * (ni + 1.);
      _an[ni] = nn;
      const double d = std::sqrt((2. * (ni + 1.) + 1.) / nn);
      _dd[ni] = d;
      for (uint_fast32_t nj = 0; nj < ni + 1; ++nj) {
        const double ddd = 0.5 * d * _dd[nj];
        _ann(ni, nj) = ddd;
        _ann(nj, ni) = ddd;
      }
    }
    SpecialFunctions::get_gauss_legendre_points_and_weigths(
        2 * ngauss, _costheta, _weights);
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      const double this_sintheta2inv =
          1. / (1. - _costheta[ig] * _costheta[ig]);
      _sintheta2inv[ig] = this_sintheta2inv;
      _sintheta2inv[2 * ngauss - ig - 1] = this_sintheta2inv;
      const double this_sinthetainv = std::sqrt(this_sintheta2inv);
      _sinthetainv[ig] = this_sinthetainv;
      _sinthetainv[2 * ngauss - ig - 1] = this_sinthetainv;
    }
    SpecialFunctions::get_r_dr_spheroid(_costheta, R_V, axis_ratio, _r2,
                                        _dr_over_r);
    const std::complex<double> mrinv = 1. / refractive_index;
    for (uint_fast32_t i = 0; i < 2 * ngauss; ++i) {
      const double r = std::sqrt(_r2[i]);
      _kr[i] = _k * r;
      _krmr[i] = refractive_index * _kr[i];
      _krinv[i] = 1. / _kr[i];
      _krmrinv[i] = mrinv * _krinv[i];
    }
    for (uint_fast32_t ig = 0; ig < 2 * ngauss; ++ig) {
      SpecialFunctions::spherical_j_jdj_array(nmax, _kr[ig], _jkr.get_row(ig),
                                              _djkr.get_row(ig));
      SpecialFunctions::spherical_y_ydy_array(nmax, _kr[ig], _ykr.get_row(ig),
                                              _dykr.get_row(ig));
      SpecialFunctions::spherical_j_jdj_array(
          nmax, _krmr[ig], _jkrmr.get_row(ig), _djkrmr.get_row(ig));
    }

    const uint_fast32_t nmax2 = 2 * nmax;

    std::vector<int_fast8_t> signs(nmax2);
    Matrix<double> wigner_d(2 * ngauss, nmax);
    Matrix<double> dwigner_d(2 * ngauss, nmax);
    std::vector<double> wr2(ngauss);
    Matrix<std::complex<double>> J12(nmax, nmax);
    Matrix<std::complex<double>> J21(nmax, nmax);
    Matrix<std::complex<double>> RgJ12(nmax, nmax);
    Matrix<std::complex<double>> RgJ21(nmax, nmax);
    Matrix<std::complex<double>> Q(nmax2, nmax2);
    Matrix<std::complex<double>> RgQ(nmax2, nmax2);

    int_fast8_t si = 1;
    for (uint_fast32_t m = 0; m < nmax2; ++m) {
      si = -si;
      signs[m] = si;
    }
    for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
      const uint_fast32_t i1 = ngauss + ig;
      const uint_fast32_t i2 = ngauss - ig + 1;
      std::vector<double> dv1(nmax), dv2(nmax);
      SpecialFunctions::wigner_dn_0m(_costheta[i1 - 1], nmax, 0, &dv1[0],
                                     &dv2[0]);
      for (uint_fast32_t n = 0; n < nmax; ++n) {
        si = signs[n];
        wigner_d(i1 - 1, n) = dv1[n];
        wigner_d(i2 - 1, n) = si * dv1[n];
        dwigner_d(i1 - 1, n) = dv2[n];
        dwigner_d(i2 - 1, n) = -si * dv2[n];
      }
    }
    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      wr2[ig] = _weights[ig] * _r2[ig];
    }
    for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
      const double an1 = _an[n1 - 1];
      for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
        const double an2 = _an[n2 - 1];

        std::complex<double> this_J12, this_J21, this_RgJ12, this_RgJ21;
        // filter out half the components because of symmetry
        if (signs[n1 + n2 - 1] > 0) {
          for (uint_fast32_t ig = 1; ig < ngauss + 1; ++ig) {
            const double wigner_n1 = wigner_d(ig - 1, n1 - 1);
            const double dwigner_n1 = dwigner_d(ig - 1, n1 - 1);
            const double wigner_n2 = wigner_d(ig - 1, n2 - 1);
            const double dwigner_n2 = dwigner_d(ig - 1, n2 - 1);

            const double a12 = wigner_n1 * dwigner_n2;
            const double a21 = dwigner_n1 * wigner_n2;
            const double a22 = dwigner_n1 * dwigner_n2;

            const double jkrn1 = _jkr(ig - 1, n1 - 1);
            const double ykrn1 = _ykr(ig - 1, n1 - 1);
            // spherical Hankel function of the first kind
            const std::complex<double> hkrn1(jkrn1, ykrn1);
            const double djkrn1 = _djkr(ig - 1, n1 - 1);
            const double dykrn1 = _dykr(ig - 1, n1 - 1);
            const std::complex<double> dhkrn1(djkrn1, dykrn1);
            const std::complex<double> jkrmrn2 = _jkrmr(ig - 1, n2 - 1);
            const std::complex<double> djkrmrn2 = _djkrmr(ig - 1, n2 - 1);

            const std::complex<double> c1 = jkrmrn2 * jkrn1;
            const std::complex<double> b1 = jkrmrn2 * hkrn1;

            const std::complex<double> c2 = jkrmrn2 * djkrn1;
            const std::complex<double> b2 = jkrmrn2 * dhkrn1;

            const double krinvi = _krinv[ig - 1];
            const std::complex<double> c3 = krinvi * c1;
            const std::complex<double> b3 = krinvi * b1;

            const std::complex<double> c4 = jkrn1 * djkrmrn2;
            const std::complex<double> b4 = hkrn1 * djkrmrn2;

            const std::complex<double> krmrinvi = _krmrinv[ig - 1];
            const std::complex<double> c5 = c1 * krmrinvi;
            const std::complex<double> b5 = b1 * krmrinvi;

            const double wr2i = wr2[ig - 1];
            const double dr_over_ri = _dr_over_r[ig - 1];

            const double f1 = wr2i * a22;
            const double f2 = wr2i * dr_over_ri * an1 * a12;
            this_J12 += f1 * b2 + f2 * b3;
            this_RgJ12 += f1 * c2 + f2 * c3;

            const double f3 = wr2i * dr_over_ri * an2 * a21;
            this_J21 += f1 * b4 + f3 * b5;
            this_RgJ21 += f1 * c4 + f3 * c5;
          }
          const double an12 = 2. * _ann(n1 - 1, n2 - 1);
          J12(n1 - 1, n2 - 1) = an12 * this_J12;
          J21(n1 - 1, n2 - 1) = an12 * this_J21;
          RgJ12(n1 - 1, n2 - 1) = an12 * this_RgJ12;
          RgJ21(n1 - 1, n2 - 1) = an12 * this_RgJ21;
        }
      }
    }
    for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
      const uint_fast32_t k1 = n1;
      const uint_fast32_t kk1 = k1 + nmax;
      for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
        const uint_fast32_t k2 = n2;
        const uint_fast32_t kk2 = k2 + nmax;

        const std::complex<double> icompl(0., 1.);
        const std::complex<double> this_J12 = -icompl * J12(n1 - 1, n2 - 1);
        const std::complex<double> this_RgJ12 = -icompl * RgJ12(n1 - 1, n2 - 1);
        const std::complex<double> this_J21 = icompl * J21(n1 - 1, n2 - 1);
        const std::complex<double> this_RgJ21 = icompl * RgJ21(n1 - 1, n2 - 1);

        Q(k1 - 1, k2 - 1) = _k2mr * this_J21 + _k2 * this_J12;
        RgQ(k1 - 1, k2 - 1) = _k2mr * this_RgJ21 + _k2 * this_RgJ12;

        Q(kk1 - 1, kk2 - 1) = _k2mr * this_J12 + _k2 * this_J21;
        RgQ(kk1 - 1, kk2 - 1) = _k2mr * this_RgJ12 + _k2 * this_RgJ21;
      }
    }

    // func_TT
    Q.plu_inverse();
    for (uint_fast32_t i = 0; i < nmax; ++i) {
      const uint_fast32_t li = (i + 1) * (i + 2) - 1;
      for (uint_fast32_t j = 0; j < nmax; ++j) {
        const uint_fast32_t lj = (j + 1) * (j + 2) - 1;
        for (uint_fast32_t k = 0; k < nmax2; ++k) {
          _T(li, lj) -= RgQ(i, k) * Q(k, j);
          _T(_Lmax + li, lj) -= RgQ(nmax + i, k) * Q(k, j);
          _T(li, _Lmax + lj) -= RgQ(i, k) * Q(k, nmax + j);
          _T(_Lmax + li, _Lmax + lj) -= RgQ(nmax + i, k) * Q(k, nmax + j);
        }
      }
    }
  }

  /**
   * @brief Compute the missing elements of the T-matrix.
   */
  inline void compute_additional_elements() {

    for (uint_fast32_t m = 1; m < _nmax; ++m) {

      const double m2 = m * m;
      const uint_fast32_t nmax2 = 2 * _nmax;
      const uint_fast32_t nm = _nmax - m + 1;
      const uint_fast32_t nm2 = 2 * nm;

      std::vector<int_fast8_t> signs(nmax2);
      Matrix<double> wigner_d(2 * _ngauss, _nmax);
      Matrix<double> dwigner_d(2 * _ngauss, _nmax);
      std::vector<double> wr2(_ngauss);
      Matrix<std::complex<double>> J11(_nmax, _nmax);
      Matrix<std::complex<double>> J12(_nmax, _nmax);
      Matrix<std::complex<double>> J21(_nmax, _nmax);
      Matrix<std::complex<double>> J22(_nmax, _nmax);
      Matrix<std::complex<double>> RgJ11(_nmax, _nmax);
      Matrix<std::complex<double>> RgJ12(_nmax, _nmax);
      Matrix<std::complex<double>> RgJ21(_nmax, _nmax);
      Matrix<std::complex<double>> RgJ22(_nmax, _nmax);
      Matrix<std::complex<double>> Q(nm2, nm2);
      Matrix<std::complex<double>> RgQ(nm2, nm2);
      std::vector<double> ds(_ngauss);
      std::vector<double> dss(_ngauss);

      int_fast8_t si = 1;
      for (uint_fast32_t m = 0; m < nmax2; ++m) {
        si = -si;
        signs[m] = si;
      }
      for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {
        const uint_fast32_t i1 = _ngauss + ig;
        const uint_fast32_t i2 = _ngauss - ig + 1;
        std::vector<double> dv1(_nmax), dv2(_nmax);
        SpecialFunctions::wigner_dn_0m(_costheta[i1 - 1], _nmax, m, &dv1[0],
                                       &dv2[0]);
        for (uint_fast32_t n = 0; n < _nmax; ++n) {
          si = signs[n];
          wigner_d(i1 - 1, n) = dv1[n];
          wigner_d(i2 - 1, n) = si * dv1[n];
          dwigner_d(i1 - 1, n) = dv2[n];
          dwigner_d(i2 - 1, n) = -si * dv2[n];
        }
      }
      for (uint_fast32_t ig = 0; ig < _ngauss; ++ig) {
        // move to class member
        wr2[ig] = _weights[ig] * _r2[ig];
        ds[ig] = _sinthetainv[ig] * m * wr2[ig];
        dss[ig] = _sintheta2inv[ig] * m2;
      }
      for (uint_fast32_t n1 = m; n1 < _nmax + 1; ++n1) {
        const double an1 = _an[n1 - 1];
        for (uint_fast32_t n2 = m; n2 < _nmax + 1; ++n2) {
          const double an2 = _an[n2 - 1];

          std::complex<double> this_J11, this_J12, this_J21, this_J22,
              this_RgJ11, this_RgJ12, this_RgJ21, this_RgJ22;
          const int_fast8_t si = signs[n1 + n2 - 1];
          for (uint_fast32_t ig = 0; ig < _ngauss; ++ig) {
            const double wigner_n1 = wigner_d(ig, n1 - 1);
            const double dwigner_n1 = dwigner_d(ig, n1 - 1);
            const double wigner_n2 = wigner_d(ig, n2 - 1);
            const double dwigner_n2 = dwigner_d(ig, n2 - 1);

            const double a11 = wigner_n1 * wigner_n2;
            const double a12 = wigner_n1 * dwigner_n2;
            const double a21 = dwigner_n1 * wigner_n2;

            const double jkrn1 = _jkr(ig, n1 - 1);
            const double ykrn1 = _ykr(ig, n1 - 1);
            // spherical Hankel function of the first kind
            const std::complex<double> hkrn1(jkrn1, ykrn1);
            const double djkrn1 = _djkr(ig, n1 - 1);
            const double dykrn1 = _dykr(ig, n1 - 1);
            const std::complex<double> dhkrn1(djkrn1, dykrn1);
            const std::complex<double> jkrmrn2 = _jkrmr(ig, n2 - 1);
            const std::complex<double> djkrmrn2 = _djkrmr(ig, n2 - 1);

            const std::complex<double> c1 = jkrmrn2 * jkrn1;
            const std::complex<double> b1 = jkrmrn2 * hkrn1;

            const std::complex<double> c2 = jkrmrn2 * djkrn1;
            const std::complex<double> b2 = jkrmrn2 * dhkrn1;

            const double krinvi = _krinv[ig];

            const std::complex<double> c4 = djkrmrn2 * jkrn1;
            const std::complex<double> b4 = djkrmrn2 * hkrn1;

            const std::complex<double> krmrinvi = _krmrinv[ig];

            const double dr_over_ri = _dr_over_r[ig];

            if (si < 0) {
              const double dsi = ds[ig];
              const double aa1 = a12 + a21;

              const std::complex<double> c6 = djkrmrn2 * djkrn1;
              const std::complex<double> b6 = djkrmrn2 * dhkrn1;

              const std::complex<double> c7 = c4 * krinvi;
              const std::complex<double> b7 = b4 * krinvi;

              const std::complex<double> c8 = c2 * krmrinvi;
              const std::complex<double> b8 = b2 * krmrinvi;
              const double e1 = dsi * aa1;
              this_J11 += e1 * b1;
              this_RgJ11 += e1 * c1;

              const double e2 = dsi * dr_over_ri * a11 * an1;
              const double e3 = dsi * dr_over_ri * a11 * an2;
              this_J22 += e1 * b6 + e2 * b7 + e3 * b8;
              this_RgJ22 += e1 * c6 + e2 * c7 + e3 * c8;
            } else {
              const double wr2i = wr2[ig];
              const double a22 = dwigner_n1 * dwigner_n2;
              const double aa2 = a11 * dss[ig] + a22;
              const std::complex<double> c3 = krinvi * c1;
              const std::complex<double> b3 = krinvi * b1;
              const std::complex<double> c5 = c1 * krmrinvi;
              const std::complex<double> b5 = b1 * krmrinvi;

              const double f1 = wr2i * aa2;
              const double f2 = wr2i * dr_over_ri * an1 * a12;
              this_J12 += f1 * b2 + f2 * b3;
              this_RgJ12 += f1 * c2 + f2 * c3;

              const double f3 = wr2i * dr_over_ri * an2 * a21;
              this_J21 += f1 * b4 + f3 * b5;
              this_RgJ21 += f1 * c4 + f3 * c5;
            }
          }
          const double an12 = 2 * _ann(n1 - 1, n2 - 1);
          J11(n1 - 1, n2 - 1) = this_J11 * an12;
          J12(n1 - 1, n2 - 1) = this_J12 * an12;
          J21(n1 - 1, n2 - 1) = this_J21 * an12;
          J22(n1 - 1, n2 - 1) = this_J22 * an12;
          RgJ11(n1 - 1, n2 - 1) = this_RgJ11 * an12;
          RgJ12(n1 - 1, n2 - 1) = this_RgJ12 * an12;
          RgJ21(n1 - 1, n2 - 1) = this_RgJ21 * an12;
          RgJ22(n1 - 1, n2 - 1) = this_RgJ22 * an12;
        }
      }
      for (uint_fast32_t n1 = m; n1 < _nmax + 1; ++n1) {
        const uint_fast32_t k1 = n1 - m + 1;
        const uint_fast32_t kk1 = k1 + nm;
        for (uint_fast32_t n2 = m; n2 < _nmax + 1; ++n2) {
          const uint_fast32_t k2 = n2 - m + 1;
          const uint_fast32_t kk2 = k2 + nm;

          const std::complex<double> icompl(0., 1.);
          const std::complex<double> this_J11 = -J11(n1 - 1, n2 - 1);
          const std::complex<double> this_RgJ11 = -RgJ11(n1 - 1, n2 - 1);
          const std::complex<double> this_J12 = -icompl * J12(n1 - 1, n2 - 1);
          const std::complex<double> this_RgJ12 =
              -icompl * RgJ12(n1 - 1, n2 - 1);
          const std::complex<double> this_J21 = icompl * J21(n1 - 1, n2 - 1);
          const std::complex<double> this_RgJ21 =
              icompl * RgJ21(n1 - 1, n2 - 1);
          const std::complex<double> this_J22 = -J22(n1 - 1, n2 - 1);
          const std::complex<double> this_RgJ22 = -RgJ22(n1 - 1, n2 - 1);

          Q(k1 - 1, k2 - 1) = _k2mr * this_J21 + _k2 * this_J12;
          RgQ(k1 - 1, k2 - 1) = _k2mr * this_RgJ21 + _k2 * this_RgJ12;

          Q(k1 - 1, kk2 - 1) = _k2mr * this_J11 + _k2 * this_J22;
          RgQ(k1 - 1, kk2 - 1) = _k2mr * this_RgJ11 + _k2 * this_RgJ22;

          Q(kk1 - 1, k2 - 1) = _k2mr * this_J22 + _k2 * this_J11;
          RgQ(kk1 - 1, k2 - 1) = _k2mr * this_RgJ22 + _k2 * this_RgJ11;

          Q(kk1 - 1, kk2 - 1) = _k2mr * this_J12 + _k2 * this_J21;
          RgQ(kk1 - 1, kk2 - 1) = _k2mr * this_RgJ12 + _k2 * this_RgJ21;
        }
      }

      Q.plu_inverse();
      for (uint_fast32_t i = 0; i < nm; ++i) {
        const uint_fast32_t lip = (m + i) * (m + i + 1) + m - 1;
        const uint_fast32_t lim = (m + i) * (m + i + 1) - m - 1;
        for (uint_fast32_t j = 0; j < nm; ++j) {
          const uint_fast32_t ljp = (m + j) * (m + j + 1) + m - 1;
          const uint_fast32_t ljm = (m + j) * (m + j + 1) - m - 1;
          for (uint_fast32_t k = 0; k < nm2; ++k) {
            _T(lip, ljp) -= RgQ(i, k) * Q(k, j);
            _T(_Lmax + lip, ljp) -= RgQ(nm + i, k) * Q(k, j);
            _T(lip, _Lmax + ljp) -= RgQ(i, k) * Q(k, nm + j);
            _T(_Lmax + lip, _Lmax + ljp) -= RgQ(nm + i, k) * Q(k, nm + j);

            _T(lim, ljm) -= RgQ(i, k) * Q(k, j);
            _T(_Lmax + lim, ljm) -= RgQ(nm + i, k) * Q(k, j);
            _T(lim, _Lmax + ljm) -= RgQ(i, k) * Q(k, nm + j);
            _T(_Lmax + lim, _Lmax + ljm) -= RgQ(nm + i, k) * Q(k, nm + j);
          }
        }
      }
    }
  }

  /**
   * @brief Get a specific element of the T-matrix.
   *
   * The T-matrix consists of four blocks:
   * @f[
   *    T = \begin{pmatrix}
   *      T^{(11)} & T^{(12)} \\
   *      T^{(21)} & T^{(22)}
   *    \end{pmatrix}.
   * @f]
   * Each of these blocks is an @f$L_{max}\times{}L_{max}@f$ matrix indexed
   * using the combined index
   * @f[
   *    l = n (n+1) + m,
   * @f]
   * where @f$n = 1...n_{max}@f$ and @f$m=-n...n@f$. This notation is based on
   * Tsang, Kong & Shin, 1984, Radio Science, 19, 629
   * (https://doi.org/10.1029/RS019i002p00629).
   *
   * This function returns the element @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   *
   * @param i1 Row index of the desired T-matrix quarter, @f$i_1@f$.
   * @param n1 Order of the row index of the element, @f$n_1@f$.
   * @param m1 Degree of the row index of the element, @f$m_1@f$.
   * @param i2 Column index of the desired T-matrix quarter, @f$i_2@f$.
   * @param n2 Order of the column index of the element, @f$n_2@f$.
   * @param m2 Degree of the column index of the element, @f$m_2@f$.
   * @return Corresponding element, @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   */
  inline const std::complex<double> &
  operator()(const uint_fast8_t i1, const uint_fast32_t n1,
             const int_fast32_t m1, const uint_fast8_t i2,
             const uint_fast32_t n2, const int_fast32_t m2) const {

    const uint_fast32_t l1 = i1 * _Lmax + n1 * (n1 + 1) + m1 - 1;
    const uint_fast32_t l2 = i2 * _Lmax + n2 * (n2 + 1) + m2 - 1;
    return _T(l1, l2);
  }
};

#endif // TMATRIX_HPP
