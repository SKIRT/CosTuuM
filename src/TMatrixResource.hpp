/**
 * @file TMatrixResource.hpp
 *
 * @brief TMatrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TMATRIXRESOURCE_HPP
#define TMATRIXRESOURCE_HPP

#include "Configuration.hpp"
#include "GaussBasedResources.hpp"
#include "InteractionResource.hpp"
#include "Matrix.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "WignerDResources.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/// predeclare tasks
class TMatrixM0Task;
class TMatrixMallTask;

/**
 * @brief TMatrix.
 */
class TMatrixResource : public Resource {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

  /*! @brief Grant access to @f$m>0@f$ computation task. */
  friend class TMatrixMallTask;

private:
  /*! @brief Actual size of the matrix stored in this resource. */
  uint_fast32_t _nmax;

  /*! @brief T-matrix itself. Is in fact a @f$n_{max}+1@f$ element vector for
   *  which every element is a @f$2n_{max}\times{}2n_{max}@f$ matrix. */
  std::vector<Matrix<std::complex<float_type>>> _T;

public:
  /**
   * @brief Constructor.
   *
   * @param maximum_order Maximum order of T matrix that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline TMatrixResource(const uint_fast32_t maximum_order) : _nmax(0) {
    _T.reserve(maximum_order + 1);
    for (uint_fast32_t m = 0; m < maximum_order + 1; ++m) {
      const uint_fast32_t nm = maximum_order + 1 - m;
      _T.push_back(Matrix<std::complex<float_type>>(2 * nm, 2 * nm));
    }
  }
};

/**
 * @brief Auxiliary data space for the T matrix calculation.
 */
class TMatrixAuxiliarySpace : public Resource {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

private:
  /*! @brief Q factor terms involving MM interactions. */
  Matrix<std::complex<float_type>> _J11;

  /*! @brief Q factor terms involving MN interactions. */
  Matrix<std::complex<float_type>> _J12;

  /*! @brief Q factor terms involving NM interactions. */
  Matrix<std::complex<float_type>> _J21;

  /*! @brief Q factor terms involving NN interactions. */
  Matrix<std::complex<float_type>> _J22;

  /*! @brief Regular Q factor terms involving MM interactions. */
  Matrix<std::complex<float_type>> _RgJ11;

  /*! @brief Regular Q factor terms involving MN interactions. */
  Matrix<std::complex<float_type>> _RgJ12;

  /*! @brief Regular Q factor terms involving NM interactions. */
  Matrix<std::complex<float_type>> _RgJ21;

  /*! @brief Regular Q factor terms involving NN interactions. */
  Matrix<std::complex<float_type>> _RgJ22;

  /*! @brief Q matrix. */
  Matrix<std::complex<float_type>> _Q;

  /*! @brief Regular Q matrix. */
  Matrix<std::complex<float_type>> _RgQ;

public:
  /**
   * @brief Constructor.
   *
   * @param maximum_order Maximum order of T matrix that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline TMatrixAuxiliarySpace(const uint_fast32_t maximum_order)
      : _J11(maximum_order, maximum_order), _J12(maximum_order, maximum_order),
        _J21(maximum_order, maximum_order), _J22(maximum_order, maximum_order),
        _RgJ11(maximum_order, maximum_order),
        _RgJ12(maximum_order, maximum_order),
        _RgJ21(maximum_order, maximum_order),
        _RgJ22(maximum_order, maximum_order),
        _Q(2 * maximum_order, 2 * maximum_order),
        _RgQ(2 * maximum_order, 2 * maximum_order) {}
};

/**
 * @brief Compute the @f$m=0@f$ part of the T-matrix for a given @f$n_{max}@f$.
 */
class TMatrixM0Task : public Task {
private:
  /*! @brief Maximum order of the matrix, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Precomputed @f$n@f$ factors (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Precomputed Gauss-Legendre quadrature points (read only). */
  const GaussBasedResources &_quadrature_points;

  /*! @brief Precomputed geometry specific quadrature points (read only). */
  const ParticleGeometryResource &_geometry;

  /*! @brief Precomputed interaction specific quadrature points (read only). */
  const InteractionResource &_interaction;

  /*! @brief Precomputed @f$m=0@f$ Wigner D functions (read only). */
  const WignerDResources &_wigner;

  /*! @brief Auxiliary space used to store intermediate calculations. */
  TMatrixAuxiliarySpace &_aux;

  /*! @brief TMatrix space containing the result. */
  TMatrixResource &_Tmatrix;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order of the matrix, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param nfactors Precomputed @f$n@f$ factors (read only).
   * @param quadrature_points Precomputed Gauss-Legendre quadrature points (read
   * only).
   * @param geometry Precomputed geometry specific quadrature points (read
   * only).
   * @param interaction Precomputed interaction specific quadrature points (read
   * only).
   * @param wigner Precomputed @f$m=0@f$ Wigner D functions (read only).
   * @param aux Auxiliary space used to store intermediate calculations.
   * @param Tmatrix TMatrix space containing the result.
   */
  inline TMatrixM0Task(const uint_fast32_t nmax, const uint_fast32_t ngauss,
                       const NBasedResources &nfactors,
                       const GaussBasedResources &quadrature_points,
                       const ParticleGeometryResource &geometry,
                       const InteractionResource &interaction,
                       const WignerDResources &wigner,
                       TMatrixAuxiliarySpace &aux, TMatrixResource &Tmatrix)
      : _nmax(nmax), _ngauss(ngauss), _nfactors(nfactors),
        _quadrature_points(quadrature_points), _geometry(geometry),
        _interaction(interaction), _wigner(wigner), _aux(aux),
        _Tmatrix(Tmatrix) {}

  /**
   * @brief Compute the @f$m=0@f$ elements of the T-matrix.
   */
  virtual void execute() {
    int_fast8_t sign_outer = 1;
    for (uint_fast32_t n1 = 1; n1 < _nmax + 1; ++n1) {
      sign_outer = -sign_outer;
      int_fast8_t sign_inner = sign_outer;
      // n1 * (n1 + 1)
      const float_type n1n1p1 = n1 * (n1 + 1.);
      for (uint_fast32_t n2 = 1; n2 < _nmax + 1; ++n2) {
        sign_inner = -sign_inner;
        // n2 * (n2 + 1)
        const float_type n2n2p1 = n2 * (n2 + 1.);

        std::complex<float_type> this_J12, this_J21, this_RgJ12, this_RgJ21;
        // filter out half the components because of symmetry
        if (sign_inner > 0) {
          for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {
            const float_type wigner_n1 = _wigner.get_wigner_d(ig - 1, n1);
            const float_type dwigner_n1 = _wigner.get_dwigner_d(ig - 1, n1);
            const float_type wigner_n2 = _wigner.get_wigner_d(ig - 1, n2);
            const float_type dwigner_n2 = _wigner.get_dwigner_d(ig - 1, n2);

            const float_type wn1dwn2 = wigner_n1 * dwigner_n2;
            const float_type dwn1wn2 = dwigner_n1 * wigner_n2;
            const float_type dwn1dwn2 = dwigner_n1 * dwigner_n2;

            const float_type jkrn1 = _interaction.get_jkr(ig - 1, n1 - 1);
            const float_type ykrn1 = _interaction.get_ykr(ig - 1, n1 - 1);
            // spherical Hankel function of the first kind
            const std::complex<float_type> hkrn1(jkrn1, ykrn1);
            const float_type djkrn1 = _interaction.get_djkr(ig - 1, n1 - 1);
            const float_type dykrn1 = _interaction.get_dykr(ig - 1, n1 - 1);
            const std::complex<float_type> dhkrn1(djkrn1, dykrn1);
            const std::complex<float_type> jkrmrn2 =
                _interaction.get_jkrmr(ig - 1, n2 - 1);
            const std::complex<float_type> djkrmrn2 =
                _interaction.get_djkrmr(ig - 1, n2 - 1);

            const std::complex<float_type> c1 = jkrmrn2 * jkrn1;
            const std::complex<float_type> b1 = jkrmrn2 * hkrn1;

            const std::complex<float_type> c2 = jkrmrn2 * djkrn1;
            const std::complex<float_type> b2 = jkrmrn2 * dhkrn1;

            const float_type krinvi = _interaction.get_krinv(ig - 1);
            const std::complex<float_type> c3 = krinvi * c1;
            const std::complex<float_type> b3 = krinvi * b1;

            const std::complex<float_type> c4 = jkrn1 * djkrmrn2;
            const std::complex<float_type> b4 = hkrn1 * djkrmrn2;

            const std::complex<float_type> krmrinvi =
                _interaction.get_krmrinv(ig - 1);
            const std::complex<float_type> c5 = c1 * krmrinvi;
            const std::complex<float_type> b5 = b1 * krmrinvi;

            const float_type wr2i = _quadrature_points.get_weight(ig - 1) *
                                    _geometry.get_r2(ig - 1);
            const float_type dr_over_ri = _geometry.get_dr_over_r(ig - 1);

            const float_type f1 = wr2i * dwn1dwn2;
            const float_type f2 = wr2i * dr_over_ri * n1n1p1 * wn1dwn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*hn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * hn1(kr) / kr
            this_J12 += f1 * b2 + f2 * b3;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*jn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * jn1(kr) / kr
            this_RgJ12 += f1 * c2 + f2 * c3;

            const float_type f3 = wr2i * dr_over_ri * n2n2p1 * dwn1wn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * hn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * hn1(kr) / krmr
            this_J21 += f1 * b4 + f3 * b5;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * jn1(kr) / krmr
            this_RgJ21 += f1 * c4 + f3 * c5;
          }
          // prefactor sqrt{(2n1+1)*(2n2+1)/[n1*(n1+1)*n2*(n2+1)]}
          const float_type an12 = 2. * _nfactors.get_ann(n1 - 1, n2 - 1);
          _aux._J12(n1 - 1, n2 - 1) = an12 * this_J12;
          _aux._J21(n1 - 1, n2 - 1) = an12 * this_J21;
          _aux._RgJ12(n1 - 1, n2 - 1) = an12 * this_RgJ12;
          _aux._RgJ21(n1 - 1, n2 - 1) = an12 * this_RgJ21;
        }
      }
    }
    for (uint_fast32_t n1 = 1; n1 < _nmax + 1; ++n1) {
      const uint_fast32_t k1 = n1;
      const uint_fast32_t kk1 = k1 + _nmax;
      for (uint_fast32_t n2 = 1; n2 < _nmax + 1; ++n2) {
        const uint_fast32_t k2 = n2;
        const uint_fast32_t kk2 = k2 + _nmax;

        const std::complex<float_type> icompl(0., 1.);
        // no idea why we multiply with i: completely unnecessary...
        // (code also works if you leave out the i factor)
        // sign differences are due to a sign difference between the
        // implementation and documentation
        const std::complex<float_type> this_J12 =
            -icompl * _aux._J12(n1 - 1, n2 - 1);
        const std::complex<float_type> this_RgJ12 =
            -icompl * _aux._RgJ12(n1 - 1, n2 - 1);
        const std::complex<float_type> this_J21 =
            icompl * _aux._J21(n1 - 1, n2 - 1);
        const std::complex<float_type> this_RgJ21 =
            icompl * _aux._RgJ21(n1 - 1, n2 - 1);

        _aux._Q(k1 - 1, k2 - 1) = _interaction.get_k2mr() * this_J21 +
                                  _interaction.get_k2() * this_J12;
        _aux._RgQ(k1 - 1, k2 - 1) = _interaction.get_k2mr() * this_RgJ21 +
                                    _interaction.get_k2() * this_RgJ12;

        _aux._Q(kk1 - 1, kk2 - 1) = _interaction.get_k2mr() * this_J12 +
                                    _interaction.get_k2() * this_J21;
        _aux._RgQ(kk1 - 1, kk2 - 1) = _interaction.get_k2mr() * this_RgJ12 +
                                      _interaction.get_k2() * this_RgJ21;
      }
    }

    // func_TT
    _aux._Q.plu_inverse();

    const uint_fast32_t nmax2 = 2 * _nmax;
    for (uint_fast32_t i = 0; i < _nmax; ++i) {
      for (uint_fast32_t j = 0; j < _nmax; ++j) {
        for (uint_fast32_t k = 0; k < nmax2; ++k) {
          _Tmatrix._T[0](i, j) -= _aux._RgQ(i, k) * _aux._Q(k, j);
          _Tmatrix._T[0](_nmax + i, j) -=
              _aux._RgQ(_nmax + i, k) * _aux._Q(k, j);
          _Tmatrix._T[0](i, _nmax + j) -=
              _aux._RgQ(i, k) * _aux._Q(k, _nmax + j);
          _Tmatrix._T[0](_nmax + i, _nmax + j) -=
              _aux._RgQ(_nmax + i, k) * _aux._Q(k, _nmax + j);
        }
      }
    }
  }
};

#endif // TMATRIXRESOURCE_HPP
