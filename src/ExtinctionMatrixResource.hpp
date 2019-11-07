/**
 * @file ExtinctionMatrixResource.hpp
 *
 * @brief Extinction matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef EXTINCTIONMATRIXRESOURCE_HPP
#define EXTINCTIONMATRIXRESOURCE_HPP

#include "Configuration.hpp"
#include "Matrix.hpp"
#include "QuickSchedWrapper.hpp"
#include "TMatrixResource.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Extinction matrix.
 */
class ExtinctionMatrixResource : public Task,
                                 public Resource,
                                 public Computable {
private:
  /*! @brief Input zenith angle (in radians). */
  const float_type _theta;

  /*! @brief Input azimuth angle (in radians). */
  const float_type _phi;

  /*! @brief T-matrix to use. */
  const TMatrixResource &_Tmatrix;

  /*! @brief Extinction matrix. */
  Matrix<float_type> _K;

public:
  /**
   * @brief Constructor.
   *
   * @param theta Input zenith angle (in radians).
   * @param phi Input azimuth angle (in radians).
   * @param Tmatrix T-matrix to use.
   */
  ExtinctionMatrixResource(const float_type theta, const float_type phi,
                           const TMatrixResource &Tmatrix)
      : _theta(theta), _phi(phi), _Tmatrix(Tmatrix), _K(4, 4) {}

  virtual ~ExtinctionMatrixResource() {}

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
   * @brief Compute the extinction matrix.
   */
  virtual void execute() {

    // compute all sines and cosines in one go; we need all of them anyway
    const float_type costheta_l = cos(_theta);
    const float_type sintheta_l = sin(_theta);
    const float_type cosphi_l = cos(_phi);
    const float_type sinphi_l = sin(_phi);

    const float_type cosphirel_in = cosphi_l;
    const float_type sinphirel_in = sinphi_l;

    const float_type costheta_p_in = costheta_l;
    const float_type sintheta_p_in = sintheta_l;
    float_type cosphi_p_in, sinphi_p_in, sintheta_p_in_inv;
    if (sintheta_p_in != 0.) {
      sintheta_p_in_inv = 1. / sintheta_p_in;
      cosphi_p_in = (sintheta_l * cosphirel_in) * sintheta_p_in_inv;
      sinphi_p_in = sintheta_l * sinphirel_in * sintheta_p_in_inv;
    } else {
      cosphi_p_in = 1.;
      sinphi_p_in = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_in_inv = 9000.;
    }

    const float_type cosphirel_out = cosphi_l;
    const float_type sinphirel_out = sinphi_l;

    const float_type costheta_p_out = costheta_l;
    const float_type sintheta_p_out = sintheta_l;
    float_type cosphi_p_out, sinphi_p_out, sintheta_p_out_inv;
    if (sintheta_p_out != 0.) {
      sintheta_p_out_inv = 1. / sintheta_p_out;
      cosphi_p_out = sintheta_l * cosphirel_out * sintheta_p_out_inv;
      sinphi_p_out = sintheta_l * sinphirel_out * sintheta_p_out_inv;
    } else {
      cosphi_p_out = 1.;
      sinphi_p_out = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_out_inv = 9000.;
    }

    Matrix<float_type> AL_in(3, 2);
    AL_in(0, 0) = costheta_l * cosphi_l;
    AL_in(0, 1) = -sinphi_l;
    AL_in(1, 0) = costheta_l * sinphi_l;
    AL_in(1, 1) = cosphi_l;
    AL_in(2, 0) = -sintheta_l;
    // AL_in(2,1) remains 0.

    Matrix<float_type> AP_in(2, 3);
    AP_in(0, 0) = costheta_p_in * cosphi_p_in;
    AP_in(0, 1) = costheta_p_in * sinphi_p_in;
    AP_in(0, 2) = -sintheta_p_in;
    AP_in(1, 0) = -sinphi_p_in;
    AP_in(1, 1) = cosphi_p_in;
    // AP_in(1,2) remains 0.

    Matrix<float_type> AL_out(3, 2);
    AL_out(0, 0) = costheta_l * cosphi_l;
    AL_out(0, 1) = -sinphi_l;
    AL_out(1, 0) = costheta_l * sinphi_l;
    AL_out(1, 1) = cosphi_l;
    AL_out(2, 0) = -sintheta_l;
    // AL_out(2,1) remains 0.

    Matrix<float_type> AP_out(2, 3);
    AP_out(0, 0) = costheta_p_out * cosphi_p_out;
    AP_out(0, 1) = costheta_p_out * sinphi_p_out;
    AP_out(0, 2) = -sintheta_p_out;
    AP_out(1, 0) = -sinphi_p_out;
    AP_out(1, 1) = cosphi_p_out;
    // AP_out(1,2) remains 0.

    Matrix<float_type> R_in(2, 2);
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          R_in(i, j) += AP_in(i, k) * AL_in(k, j);
        }
      }
    }

    Matrix<float_type> R_out(2, 2);
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          R_out(i, j) += AP_out(i, k) * AL_out(k, j);
        }
      }
    }

    // manually invert the 2x2 matrix R_out
    const float_type d =
        1. / (R_out(0, 0) * R_out(1, 1) - R_out(0, 1) * R_out(1, 0));
    const float_type temp = R_out(0, 0);
    R_out(0, 0) = R_out(1, 1) * d;
    R_out(0, 1) = -R_out(0, 1) * d;
    R_out(1, 0) = -R_out(1, 0) * d;
    R_out(1, 1) = temp * d;

    // precompute the c factors
    const std::complex<float_type> icompl(0., 1.);
    Matrix<std::complex<float_type>> c(_Tmatrix.get_nmax(),
                                       _Tmatrix.get_nmax());
    std::complex<float_type> icomp_pow_nn = icompl;
    for (uint_fast32_t nn = 1; nn < _Tmatrix.get_nmax() + 1; ++nn) {
      std::complex<float_type> icomp_pow_m_n_m_1(-1.);
      for (uint_fast32_t n = 1; n < _Tmatrix.get_nmax() + 1; ++n) {
        // icomp_pow_nn*icomp_pow_m_n_m_1 now equals i^(nn - n - 1)
        c(n - 1, nn - 1) = icomp_pow_m_n_m_1 * icomp_pow_nn *
                           float_type(sqrt((2. * n + 1.) * (2. * nn + 1.) /
                                           (n * nn * (n + 1.) * (nn + 1.))));
        icomp_pow_m_n_m_1 /= icompl;
      }
      icomp_pow_nn *= icompl;
    }

    // now compute the matrix S^P
    // we precompute e^{i(phi_out-phi_in)}
    const std::complex<float_type> expiphi_p_out_m_in(
        cosphi_p_out * cosphi_p_in + sinphi_p_out * sinphi_p_in,
        sinphi_p_out * cosphi_p_in - cosphi_p_out * sinphi_p_in);
    // e^{im(phi_out-phi_in)} is computed recursively, starting with the value
    // for m=0: 1
    std::complex<float_type> expimphi_p_out_m_in(1., 0.);
    Matrix<std::complex<float_type>> S(2, 2);
    // instead of summing over n and n', we sum over m, since then we can reuse
    // the e^{im(phi_out-phi_in)}, pi and tau factors
    for (uint_fast32_t m = 0; m < _Tmatrix.get_nmax() + 1; ++m) {
      // only n and n' values larger than or equal to m have non-trivial
      // contributions to the S matrix
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));

      // precompute the pi and tau functions for this value of m
      std::vector<float_type> pi_in(_Tmatrix.get_nmax()),
          tau_in(_Tmatrix.get_nmax());
      SpecialFunctions::wigner_dn_0m_sinx(
          costheta_p_in, sintheta_p_in, sintheta_p_in_inv, _Tmatrix.get_nmax(),
          m, &pi_in[0], &tau_in[0]);
      std::vector<float_type> pi_out(_Tmatrix.get_nmax()),
          tau_out(_Tmatrix.get_nmax());
      SpecialFunctions::wigner_dn_0m_sinx(
          costheta_p_out, sintheta_p_out, sintheta_p_out_inv,
          _Tmatrix.get_nmax(), m, &pi_out[0], &tau_out[0]);

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      // now perform the actual sums over n and n'
      for (uint_fast32_t nn = nmin; nn < _Tmatrix.get_nmax() + 1; ++nn) {

        // get the specific pi and tau for this n'
        const float_type pi_nn = m * pi_in[nn - 1];
        const float_type tau_nn = tau_in[nn - 1];

        for (uint_fast32_t n = nmin; n < _Tmatrix.get_nmax() + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n = m * pi_out[n - 1];
          const float_type tau_n = tau_out[n - 1];

          // get the c factor for these values of n and n'
          const std::complex<float_type> c_nnn = c(n - 1, nn - 1);

          // get the T11 and T22 elements for this m, n and n' (we need these
          // in all cases)
          const std::complex<float_type> T11nmnnm = _Tmatrix(0, n, m, 0, nn, m);
          const std::complex<float_type> T22nmnnm = _Tmatrix(1, n, m, 1, nn, m);
          // if m=0, the T12 and T21 matrices are trivially zero, and we can
          // simplify the expression for S
          if (m == 0) {
            const std::complex<float_type> factor = c_nnn * tau_n * tau_nn;
            S(0, 0) += factor * T22nmnnm;
            S(1, 1) += factor * T11nmnnm;
          } else {
            // in the general case m=/=0, we also need the T12 and T21 elements
            // for this m, n and n'
            const std::complex<float_type> T12nmnnm =
                _Tmatrix(0, n, m, 1, nn, m);
            const std::complex<float_type> T21nmnnm =
                _Tmatrix(1, n, m, 0, nn, m);

            // due to m symmetry, S11 and S22 only have the cosine factor,
            // while S12 and S21 only have the sine factor
            const std::complex<float_type> real_factor = c_nnn * fcos;
            const std::complex<float_type> imag_factor = c_nnn * fsin;

            // precompute the pi and tau factor combinations
            const float_type pi_pi = pi_n * pi_nn;
            const float_type pi_tau = pi_n * tau_nn;
            const float_type tau_pi = tau_n * pi_nn;
            const float_type tau_tau = tau_n * tau_nn;

            S(0, 0) += real_factor * (T11nmnnm * pi_pi + T21nmnnm * tau_pi +
                                      T12nmnnm * pi_tau + T22nmnnm * tau_tau);
            S(0, 1) += imag_factor * (T11nmnnm * pi_tau + T21nmnnm * tau_tau +
                                      T12nmnnm * pi_pi + T22nmnnm * tau_pi);
            S(1, 0) -= imag_factor * (T11nmnnm * tau_pi + T21nmnnm * pi_pi +
                                      T12nmnnm * tau_tau + T22nmnnm * pi_tau);
            S(1, 1) += real_factor * (T11nmnnm * tau_tau + T21nmnnm * pi_tau +
                                      T12nmnnm * tau_pi + T22nmnnm * pi_pi);
          }
        }
      }
    }
    // now divide all expressions by the wavenumber
    const float_type kinv = 1. / _Tmatrix.get_wavenumber();
    S(0, 0) *= kinv;
    S(0, 1) *= kinv;
    S(1, 0) *= kinv;
    S(1, 1) *= kinv;

    // perform the double 2x2 matrix product to convert S^P to S^L
    const std::complex<float_type> cS11 =
        S(0, 0) * R_in(0, 0) + S(0, 1) * R_in(1, 0);
    const std::complex<float_type> cS12 =
        S(0, 0) * R_in(0, 1) + S(0, 1) * R_in(1, 1);
    const std::complex<float_type> cS21 =
        S(1, 0) * R_in(0, 0) + S(1, 1) * R_in(1, 0);
    const std::complex<float_type> cS22 =
        S(1, 0) * R_in(0, 1) + S(1, 1) * R_in(1, 1);

    S(0, 0) = R_out(0, 0) * cS11 + R_out(0, 1) * cS21;
    S(0, 1) = R_out(0, 0) * cS12 + R_out(0, 1) * cS22;
    S(1, 0) = R_out(1, 0) * cS11 + R_out(1, 1) * cS21;
    S(1, 1) = R_out(1, 0) * cS12 + R_out(1, 1) * cS22;

    const float_type prefactor = 2. * M_PI * kinv;

    _K(0, 0) = prefactor * (S(0, 0) + S(1, 1)).imag();
    _K(1, 1) = _K(0, 0);
    _K(2, 2) = _K(0, 0);
    _K(3, 3) = _K(0, 0);

    _K(0, 1) = prefactor * (S(0, 0) - S(1, 1)).imag();
    _K(1, 0) = _K(0, 1);

    _K(0, 2) = -prefactor * (S(0, 1) + S(1, 0)).imag();
    _K(2, 0) = _K(0, 2);

    _K(0, 3) = prefactor * (S(1, 0) - S(0, 1)).real();
    _K(3, 0) = _K(0, 3);

    _K(1, 2) = prefactor * (S(1, 0) - S(0, 1)).imag();
    _K(2, 1) = -_K(1, 2);

    _K(1, 3) = -prefactor * (S(0, 1) + S(1, 0)).real();
    _K(3, 1) = -_K(1, 3);

    _K(2, 3) = prefactor * (S(1, 1) - S(0, 0)).real();
    _K(3, 2) = -_K(2, 3);

    make_available();
  }

  /**
   * @brief Access the given element of the extinction matrix.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Corresponding element of the extinction matrix.
   */
  inline float_type operator()(const uint_fast8_t i,
                               const uint_fast8_t j) const {

    ctm_assert(i < 4);
    ctm_assert(j < 4);

    check_use();

    return _K(i, j);
  }
};

#endif // EXTINCTIONMATRIXRESOURCE_HPP
