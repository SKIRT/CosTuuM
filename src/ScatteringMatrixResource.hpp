/**
 * @file ScatteringMatrixResource.hpp
 *
 * @brief Scattering matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SCATTERINGMATRIXRESOURCE_HPP
#define SCATTERINGMATRIXRESOURCE_HPP

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
 * @brief Scattering matrix.
 */
class ScatteringMatrixResource : public Task,
                                 public Resource,
                                 public Computable {
private:
  /*! @brief Azimuth angle of the scattering particle (in radians). */
  const float_type _alpha;

  /*! @brief Zenith angle of the scattering particle (in radians). */
  const float_type _beta;

  /*! @brief Input zenith angle (in radians). */
  const float_type _theta_in;

  /*! @brief Input azimuth angle (in radians). */
  const float_type _phi_in;

  /*! @brief Output zenith angle (in radians). */
  const float_type _theta_out;

  /*! @brief Output azimuth angle (in radians). */
  const float_type _phi_out;

  /*! @brief T-matrix to use. */
  const TMatrixResource &_Tmatrix;

  /*! @brief Scattering matrix. */
  Matrix<float_type> _Z;

  /*! @brief Temporary B matrix. */
  Matrix<float_type> _B;

  /*! @brief Temporary AL_in matrix. */
  Matrix<float_type> _AL_in;

  /*! @brief Temporary AP_in matrix. */
  Matrix<float_type> _AP_in;

  /*! @brief Temporary AL_out matrix. */
  Matrix<float_type> _AL_out;

  /*! @brief Temporary AP_out matrix. */
  Matrix<float_type> _AP_out;

  /*! @brief Temporary C matrix. */
  Matrix<float_type> _C;

  /*! @brief Temporary R_in matrix. */
  Matrix<float_type> _R_in;

  /*! @brief Temporary R_out matrix. */
  Matrix<float_type> _R_out;

  /*! @brief Prefactor matrix. */
  Matrix<std::complex<float_type>> _c;

  /*! @brief Temporary S matrix. */
  Matrix<std::complex<float_type>> _S;

  /*! @brief Temporary Wigner D arrays. */
  std::vector<float_type> _pi_in;

  /*! @brief Temporary Wigner D arrays. */
  std::vector<float_type> _tau_in;

  /*! @brief Temporary Wigner D arrays. */
  std::vector<float_type> _pi_out;

  /*! @brief Temporary Wigner D arrays. */
  std::vector<float_type> _tau_out;

public:
  /**
   * @brief Constructor.
   *
   * @param alpha Azimuth angle of the scattering particle (in radians).
   * @param beta Zenith angle of the scattering particle (in radians).
   * @param theta_in Input zenith angle (in radians).
   * @param phi_in Input azimuth angle (in radians).
   * @param theta_out Output zenith angle (in radians).
   * @param phi_out Output azimuth angle (in radians).
   * @param Tmatrix T-matrix to use.
   * @param nmax Maximum order of the T-matrix.
   */
  ScatteringMatrixResource(const float_type alpha, const float_type beta,
                           const float_type theta_in, const float_type phi_in,
                           const float_type theta_out, const float_type phi_out,
                           const TMatrixResource &Tmatrix,
                           const uint_fast32_t nmax)
      : _alpha(alpha), _beta(beta), _theta_in(theta_in), _phi_in(phi_in),
        _theta_out(theta_out), _phi_out(phi_out), _Tmatrix(Tmatrix), _Z(4, 4),
        _B(3, 3), _AL_in(3, 2), _AP_in(2, 3), _AL_out(3, 2), _AP_out(2, 3),
        _C(3, 2), _R_in(2, 2), _R_out(2, 2), _c(nmax, nmax), _S(2, 2),
        _pi_in(nmax), _tau_in(nmax), _pi_out(nmax), _tau_out(nmax) {}

  virtual ~ScatteringMatrixResource() {}

  /**
   * @brief Get the size in memory of a hypothetical ScatteringMatrixResource
   * object with the given parameters.
   *
   * @param maximum_order Maximum order of the hypothetical object.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t maximum_order) {
    // storage space for class variables
    size_t size = sizeof(ScatteringMatrixResource);
    // float_type matrices
    size += (16 + 9 + 5 * 6 + 2 * 4) * sizeof(float_type);
    // complex matrices
    size += (maximum_order * maximum_order + 4) * 2 * sizeof(float_type);
    // Wigner D arrays
    size += 4 * maximum_order * sizeof(float_type);
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
   * @brief Compute the scattering matrix.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    const float_type cosalpha = cos(_alpha);
    const float_type sinalpha = sin(_alpha);
    const float_type cosbeta = cos(_beta);
    const float_type sinbeta = sin(_beta);
    const float_type costheta_l_in = cos(_theta_in);
    const float_type sintheta_l_in = sin(_theta_in);
    const float_type costheta_l_out = cos(_theta_out);
    const float_type sintheta_l_out = sin(_theta_out);
    const float_type cosphi_l_in = cos(_phi_in);
    const float_type sinphi_l_in = sin(_phi_in);
    const float_type cosphi_l_out = cos(_phi_out);
    const float_type sinphi_l_out = sin(_phi_out);

    const float_type cosphirel_in =
        cosphi_l_in * cosalpha + sinphi_l_in * sinalpha;
    const float_type sinphirel_in =
        sinphi_l_in * cosalpha - cosphi_l_in * sinalpha;

    const float_type costheta_p_in =
        costheta_l_in * cosbeta + sintheta_l_in * sinbeta * cosphirel_in;
    const float_type sintheta_p_in =
        sqrt((1. - costheta_p_in) * (1. + costheta_p_in));
    float_type cosphi_p_in, sinphi_p_in, sintheta_p_in_inv;
    if (sintheta_p_in != 0.) {
      sintheta_p_in_inv = 1. / sintheta_p_in;
      cosphi_p_in =
          (cosbeta * sintheta_l_in * cosphirel_in - sinbeta * costheta_l_in) *
          sintheta_p_in_inv;
      sinphi_p_in = sintheta_l_in * sinphirel_in * sintheta_p_in_inv;
    } else {
      cosphi_p_in = 1.;
      sinphi_p_in = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_in_inv = 9000.;
    }

    const float_type cosphirel_out =
        cosphi_l_out * cosalpha + sinphi_l_out * sinalpha;
    const float_type sinphirel_out =
        sinphi_l_out * cosalpha - cosphi_l_out * sinalpha;

    const float_type costheta_p_out =
        costheta_l_out * cosbeta + sintheta_l_out * sinbeta * cosphirel_out;
    const float_type sintheta_p_out =
        sqrt((1. - costheta_p_out) * (1. + costheta_p_out));
    float_type cosphi_p_out, sinphi_p_out, sintheta_p_out_inv;
    if (sintheta_p_out != 0.) {
      sintheta_p_out_inv = 1. / sintheta_p_out;
      cosphi_p_out = (cosbeta * sintheta_l_out * cosphirel_out -
                      sinbeta * costheta_l_out) *
                     sintheta_p_out_inv;
      sinphi_p_out = sintheta_l_out * sinphirel_out * sintheta_p_out_inv;
    } else {
      cosphi_p_out = 1.;
      sinphi_p_out = 0.;
      // this value will not be used, but we set it to something anyway
      sintheta_p_out_inv = 9000.;
    }

    _B(0, 0) = cosalpha * cosbeta;
    _B(0, 1) = sinalpha * cosbeta;
    _B(0, 2) = -sinbeta;
    _B(1, 0) = -sinalpha;
    _B(1, 1) = cosalpha;
    // B(1,2) remains 0.
    _B(2, 0) = cosalpha * sinbeta;
    _B(2, 1) = sinalpha * sinbeta;
    _B(2, 2) = cosbeta;

    _AL_in(0, 0) = costheta_l_in * cosphi_l_in;
    _AL_in(0, 1) = -sinphi_l_in;
    _AL_in(1, 0) = costheta_l_in * sinphi_l_in;
    _AL_in(1, 1) = cosphi_l_in;
    _AL_in(2, 0) = -sintheta_l_in;
    // AL_in(2,1) remains 0.

    _AP_in(0, 0) = costheta_p_in * cosphi_p_in;
    _AP_in(0, 1) = costheta_p_in * sinphi_p_in;
    _AP_in(0, 2) = -sintheta_p_in;
    _AP_in(1, 0) = -sinphi_p_in;
    _AP_in(1, 1) = cosphi_p_in;
    // AP_in(1,2) remains 0.

    _AL_out(0, 0) = costheta_l_out * cosphi_l_out;
    _AL_out(0, 1) = -sinphi_l_out;
    _AL_out(1, 0) = costheta_l_out * sinphi_l_out;
    _AL_out(1, 1) = cosphi_l_out;
    _AL_out(2, 0) = -sintheta_l_out;
    // AL_out(2,1) remains 0.

    _AP_out(0, 0) = costheta_p_out * cosphi_p_out;
    _AP_out(0, 1) = costheta_p_out * sinphi_p_out;
    _AP_out(0, 2) = -sintheta_p_out;
    _AP_out(1, 0) = -sinphi_p_out;
    _AP_out(1, 1) = cosphi_p_out;
    // AP_out(1,2) remains 0.

    // C is a temporary matrix that contains B x AL_in
    for (uint_fast8_t i = 0; i < 3; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          _C(i, j) += _B(i, k) * _AL_in(k, j);
        }
      }
    }
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          _R_in(i, j) += _AP_in(i, k) * _C(k, j);
        }
      }
    }

    // now C will contain B x AL_out
    for (uint_fast8_t i = 0; i < 3; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        _C(i, j) = 0.;
        for (uint_fast8_t k = 0; k < 3; ++k) {
          _C(i, j) += _B(i, k) * _AL_out(k, j);
        }
      }
    }
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 2; ++j) {
        for (uint_fast8_t k = 0; k < 3; ++k) {
          _R_out(i, j) += _AP_out(i, k) * _C(k, j);
        }
      }
    }

    // manually invert the 2x2 matrix R_out
    const float_type d =
        1. / (_R_out(0, 0) * _R_out(1, 1) - _R_out(0, 1) * _R_out(1, 0));
    const float_type temp = _R_out(0, 0);
    _R_out(0, 0) = _R_out(1, 1) * d;
    _R_out(0, 1) = -_R_out(0, 1) * d;
    _R_out(1, 0) = -_R_out(1, 0) * d;
    _R_out(1, 1) = temp * d;

    // precompute the c factors
    const std::complex<float_type> icompl(0., 1.);
    std::complex<float_type> icomp_pow_nn = icompl;
    for (uint_fast32_t nn = 1; nn < _Tmatrix.get_nmax() + 1; ++nn) {
      std::complex<float_type> icomp_pow_m_n_m_1(-1.);
      for (uint_fast32_t n = 1; n < _Tmatrix.get_nmax() + 1; ++n) {
        // icomp_pow_nn*icomp_pow_m_n_m_1 now equals i^(nn - n - 1)
        _c(n - 1, nn - 1) = icomp_pow_m_n_m_1 * icomp_pow_nn *
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
    // instead of summing over n and n', we sum over m, since then we can reuse
    // the e^{im(phi_out-phi_in)}, pi and tau factors
    for (uint_fast32_t m = 0; m < _Tmatrix.get_nmax() + 1; ++m) {
      // only n and n' values larger than or equal to m have non-trivial
      // contributions to the S matrix
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));

      // precompute the pi and tau functions for this value of m
      SpecialFunctions::wigner_dn_0m_sinx(
          costheta_p_in, sintheta_p_in, sintheta_p_in_inv, _Tmatrix.get_nmax(),
          m, &_pi_in[0], &_tau_in[0]);
      SpecialFunctions::wigner_dn_0m_sinx(
          costheta_p_out, sintheta_p_out, sintheta_p_out_inv,
          _Tmatrix.get_nmax(), m, &_pi_out[0], &_tau_out[0]);

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      // now perform the actual sums over n and n'
      for (uint_fast32_t nn = nmin; nn < _Tmatrix.get_nmax() + 1; ++nn) {

        // get the specific pi and tau for this n'
        const float_type pi_nn = m * _pi_in[nn - 1];
        const float_type tau_nn = _tau_in[nn - 1];

        for (uint_fast32_t n = nmin; n < _Tmatrix.get_nmax() + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n = m * _pi_out[n - 1];
          const float_type tau_n = _tau_out[n - 1];

          // get the c factor for these values of n and n'
          const std::complex<float_type> c_nnn = _c(n - 1, nn - 1);

          // get the T11 and T22 elements for this m, n and n' (we need these
          // in all cases)
          const std::complex<float_type> T11nmnnm = _Tmatrix(0, n, m, 0, nn, m);
          const std::complex<float_type> T22nmnnm = _Tmatrix(1, n, m, 1, nn, m);
          // if m=0, the T12 and T21 matrices are trivially zero, and we can
          // simplify the expression for S
          if (m == 0) {
            const std::complex<float_type> factor = c_nnn * tau_n * tau_nn;
            _S(0, 0) += factor * T22nmnnm;
            _S(1, 1) += factor * T11nmnnm;
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

            _S(0, 0) += real_factor * (T11nmnnm * pi_pi + T21nmnnm * tau_pi +
                                       T12nmnnm * pi_tau + T22nmnnm * tau_tau);
            _S(0, 1) += imag_factor * (T11nmnnm * pi_tau + T21nmnnm * tau_tau +
                                       T12nmnnm * pi_pi + T22nmnnm * tau_pi);
            _S(1, 0) -= imag_factor * (T11nmnnm * tau_pi + T21nmnnm * pi_pi +
                                       T12nmnnm * tau_tau + T22nmnnm * pi_tau);
            _S(1, 1) += real_factor * (T11nmnnm * tau_tau + T21nmnnm * pi_tau +
                                       T12nmnnm * tau_pi + T22nmnnm * pi_pi);
          }
        }
      }
    }
    // now divide all expressions by the wavenumber
    const float_type kinv = 1. / _Tmatrix.get_wavenumber();
    _S(0, 0) *= kinv;
    _S(0, 1) *= kinv;
    _S(1, 0) *= kinv;
    _S(1, 1) *= kinv;

    // perform the double 2x2 matrix product to convert S^P to S^L
    const std::complex<float_type> cS11 =
        _S(0, 0) * _R_in(0, 0) + _S(0, 1) * _R_in(1, 0);
    const std::complex<float_type> cS12 =
        _S(0, 0) * _R_in(0, 1) + _S(0, 1) * _R_in(1, 1);
    const std::complex<float_type> cS21 =
        _S(1, 0) * _R_in(0, 0) + _S(1, 1) * _R_in(1, 0);
    const std::complex<float_type> cS22 =
        _S(1, 0) * _R_in(0, 1) + _S(1, 1) * _R_in(1, 1);

    _S(0, 0) = _R_out(0, 0) * cS11 + _R_out(0, 1) * cS21;
    _S(0, 1) = _R_out(0, 0) * cS12 + _R_out(0, 1) * cS22;
    _S(1, 0) = _R_out(1, 0) * cS11 + _R_out(1, 1) * cS21;
    _S(1, 1) = _R_out(1, 0) * cS12 + _R_out(1, 1) * cS22;

    const float_type half(0.5);
    _Z(0, 0) = (half * (_S(0, 0) * conj(_S(0, 0)) + _S(0, 1) * conj(_S(0, 1)) +
                        _S(1, 0) * conj(_S(1, 0)) + _S(1, 1) * conj(_S(1, 1))))
                   .real();
    _Z(0, 1) = (half * (_S(0, 0) * conj(_S(0, 0)) - _S(0, 1) * conj(_S(0, 1)) +
                        _S(1, 0) * conj(_S(1, 0)) - _S(1, 1) * conj(_S(1, 1))))
                   .real();
    _Z(0, 2) = (-_S(0, 0) * conj(_S(0, 1)) - _S(1, 1) * conj(_S(1, 0))).real();
    _Z(0, 3) =
        (icompl * (_S(0, 0) * conj(_S(0, 1)) - _S(1, 1) * conj(_S(1, 0))))
            .real();

    _Z(1, 0) = (half * (_S(0, 0) * conj(_S(0, 0)) + _S(0, 1) * conj(_S(0, 1)) -
                        _S(1, 0) * conj(_S(1, 0)) - _S(1, 1) * conj(_S(1, 1))))
                   .real();
    _Z(1, 1) = (half * (_S(0, 0) * conj(_S(0, 0)) - _S(0, 1) * conj(_S(0, 1)) -
                        _S(1, 0) * conj(_S(1, 0)) + _S(1, 1) * conj(_S(1, 1))))
                   .real();
    _Z(1, 2) = (-_S(0, 0) * conj(_S(0, 1)) + _S(1, 1) * conj(_S(1, 0))).real();
    _Z(1, 3) =
        (icompl * (_S(0, 0) * conj(_S(0, 1)) + _S(1, 1) * conj(_S(1, 0))))
            .real();

    _Z(2, 0) = (-_S(0, 0) * conj(_S(1, 0)) - _S(1, 1) * conj(_S(0, 1))).real();
    _Z(2, 1) = (-_S(0, 0) * conj(_S(1, 0)) + _S(1, 1) * conj(_S(0, 1))).real();
    _Z(2, 2) = (_S(0, 0) * conj(_S(1, 1)) + _S(0, 1) * conj(_S(1, 0))).real();
    _Z(2, 3) =
        (-icompl * (_S(0, 0) * conj(_S(1, 1)) + _S(1, 0) * conj(_S(0, 1))))
            .real();

    _Z(3, 0) =
        (icompl * (_S(1, 0) * conj(_S(0, 0)) + _S(1, 1) * conj(_S(0, 1))))
            .real();
    _Z(3, 1) =
        (icompl * (_S(1, 0) * conj(_S(0, 0)) - _S(1, 1) * conj(_S(0, 1))))
            .real();
    _Z(3, 2) =
        (-icompl * (_S(1, 1) * conj(_S(0, 0)) - _S(0, 1) * conj(_S(1, 0))))
            .real();
    _Z(3, 3) = (_S(1, 1) * conj(_S(0, 0)) - _S(0, 1) * conj(_S(1, 0))).real();

    make_available();
  }

  /**
   * @brief Access the given element of the scattering matrix.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Corresponding element of the scattering matrix.
   */
  inline float_type operator()(const uint_fast8_t i,
                               const uint_fast8_t j) const {

    ctm_assert(i < 4);
    ctm_assert(j < 4);

    check_use();

    return _Z(i, j);
  }
};

#endif // SCATTERINGMATRIXRESOURCE_HPP
