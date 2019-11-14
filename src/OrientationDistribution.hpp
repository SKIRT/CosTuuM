/**
 * @file OrientationDistribution.hpp
 *
 * @brief Class that contains the coefficients of the Legendre function
 * expansion of an orientation distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ORIENTATIONDISTRIBUTION_HPP
#define ORIENTATIONDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "SpecialFunctions.hpp"

/**
 * @brief Class that contains the coefficients of the Legendre function
 * expansion of an orientation distribution.
 */
class OrientationDistribution {
private:
  /*! @brief Normalisation factor for the distribution function. */
  float_type _normalisation_factor;

public:
  /**
   * @brief Virtual function that contains an actual expression for the
   * orientation distribution function as a function of the zenith angle
   * @f$\beta{}@f$.
   *
   * The zenith angle @f$\beta{}@f$ represents the angle between the average
   * orientation of the rotation axis of the spheroidal dust particles and the
   * specific orientation for which we want to know the probability.
   * A distribution function @f$p(\beta{}) = \delta{}(\beta{})@f$ e.g. would
   * correspond to perfect alignment of all dust particles.
   *
   * The additional cosine and sine arguments are provided because they often
   * appear in more complicated expressions for @f$p(\beta{})@f$, and also
   * because they are computed in the integrand anyway.
   *
   * @param beta Zenith angle @f$\beta{} \in{} [0, \pi{}]@f$.
   * @param cosbeta Cosine of the zenith angle, @f$\cos(\beta{})@f$.
   * @param sinbeta Sine of the zenith angle, @f$\sin(\beta{})@f$.
   * @return Orientation distribution function for that value of
   * @f$\beta{}@f$, @f$p(\beta{})@f$.
   */
  virtual float_type operator()(const float_type beta, const float_type cosbeta,
                                const float_type sinbeta) const {

    return 0.5;
  }

  /**
   * @brief Wrapper around the distribution function that does not require the
   * sine and cosine of the zenith angle.
   *
   * Provided in case the distribution function is used outside the integrand.
   *
   * @param beta Zenith angle, @f$\beta{} \in{} [0, \pi{}]@f$.
   * @return Orientation distribution function for that value of
   * @f$\beta{}@f$, @f$p(\beta{})@f$.
   */
  float_type operator()(const float_type beta) const {

    const float_type cosbeta = cos(beta);
    const float_type sinbeta = sin(beta);
    return _normalisation_factor * (*this)(beta, cosbeta, sinbeta);
  }

private:
  /*! @brief Coefficients. */
  std::vector<float_type> _coefficients;

  /**
   * @brief Auxiliary class used to wrap integrand function arguments.
   */
  class IntegrandArguments {
  private:
    /*! @brief Order of the Wigner D function in the integrand. */
    const uint_fast32_t _order;

    /*! @brief Reference to the actual OrientationDistribution object. */
    const OrientationDistribution &_orientation_distribution;

  public:
    /**
     * @brief Constructor.
     *
     * @param order Order of the Wigner D function in the integrand.
     * @param orientation_distribution Reference to the orientation
     * distribution that needs to be used.
     */
    inline IntegrandArguments(
        const uint_fast32_t order,
        const OrientationDistribution &orientation_distribution)
        : _order(order), _orientation_distribution(orientation_distribution) {}

    /**
     * @brief Get the order of the Wigner D function in the integrand.
     *
     * @return Order of the Wigner D function in the integrand.
     */
    inline uint_fast32_t get_order() const { return _order; }

    /**
     * @brief Get the orientation distribution that is used in the integrand.
     *
     * @return OrientationDistribution that is used.
     */
    inline const OrientationDistribution &get_orientation_distribution() const {
      return _orientation_distribution;
    }
  };

  /**
   * @brief Integrand for coefficients.
   *
   * @param beta Angle, @f$\beta{}@f$.
   * @param args Additional arguments: the order of the Wigner D function in the
   * integrand and a reference to the orientation distribution itself.
   * @return Value of the integrand.
   */
  static float_type integrand(const float_type beta, void *args) {

    const IntegrandArguments &arguments =
        *reinterpret_cast<IntegrandArguments *>(args);
    const uint_fast32_t order = arguments.get_order();
    const OrientationDistribution &orientation_distribution =
        arguments.get_orientation_distribution();

    const float_type cosbeta = cos(beta);
    const float_type sinbeta = sin(beta);
    float_type wigner_d_n;
    if (order > 0) {
      std::vector<float_type> wigner_d(order), dwigner_d(order);
      SpecialFunctions::wigner_dn_0m(cosbeta, order, 0, &wigner_d[0],
                                     &dwigner_d[0]);
      wigner_d_n = wigner_d[order - 1];
    } else {
      wigner_d_n = 1.;
    }
    const float_type pbeta = orientation_distribution._normalisation_factor *
                             orientation_distribution(beta, cosbeta, sinbeta);
    return sinbeta * pbeta * wigner_d_n;
  }

  /**
   * @brief Integrand used in the normalisation integral for @f$p(\beta{})@f$.
   *
   * @param beta Angle @f$\beta{}@f$.
   * @param args Additional arguments: reference to the orientation
   * distribution itself.
   * @return Value of the integrand, @f$\sin(\beta{})p(\beta{})@f$.
   */
  static float_type normalisation_integrand(const float_type beta, void *args) {
    const IntegrandArguments &arguments =
        *reinterpret_cast<IntegrandArguments *>(args);
    const OrientationDistribution &orientation_distribution =
        arguments.get_orientation_distribution();

    const float_type cosbeta = cos(beta);
    const float_type sinbeta = sin(beta);
    const float_type pbeta = orientation_distribution(beta, cosbeta, sinbeta);
    return sinbeta * pbeta;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   */
  inline OrientationDistribution(const uint_fast32_t nmax)
      : _normalisation_factor(1.), _coefficients(nmax + 1, 0.) {}

  /**
   * @brief Initialise the coefficients of the distribution.
   */
  inline void initialise() {

    // compute the normalisation factor
    {
      IntegrandArguments arguments(0, *this);
      const float_type norm =
          SpecialFunctions::gauss_legendre_quadrature<float_type>(
              normalisation_integrand, 0., M_PI, &arguments, 100, 1000, 1.e-10,
              1.e-5);
      _normalisation_factor = 1. / norm;
    }

    for (uint_fast32_t i = 0; i < _coefficients.size(); ++i) {
      IntegrandArguments arguments(i, *this);
      _coefficients[i] =
          SpecialFunctions::gauss_legendre_quadrature<float_type>(
              integrand, 0., M_PI, &arguments, 2 * (i + 1), 100 * (i + 1),
              1.e-10, 1.e-5);
    }
  }

  /**
   * @brief Get coefficient @f$p_n@f$.
   *
   * @param n Order @f$n@f$.
   * @return Corresponding coefficient.
   */
  inline float_type get_coefficient(const uint_fast32_t n) const {
    return _coefficients[n];
  }

  /**
   * @brief Get the highest order of coefficients stored in this object.
   *
   * @return Maximum order that is stored in this object.
   */
  inline uint_fast32_t get_maximum_order() const {
    return _coefficients.size() - 1;
  }
};

#endif // ORIENTATIONDISTRIBUTION_HPP
