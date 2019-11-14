/**
 * @file MishchenkoOrientationDistribution.hpp
 *
 * @brief Mishchenko (1991) imperfect alignment orientation distribution.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef MISHCHENKOORIENTATIONDISTRIBUTION_HPP
#define MISHCHENKOORIENTATIONDISTRIBUTION_HPP

#include "OrientationDistribution.hpp"

/**
 * @brief Mishchenko (1991) imperfect alignment orientation distribution.
 *
 * Based on Mishchenko, 1991, ApJ, 367, 561 (https://doi.org/10.1086/169652).
 *
 * The alignment distribution is given by
 * @f[
 *  p(\beta{}) = C\left[
 *    p_1 + \frac{5}{4} p_2 \left(\cos^2(\beta{}) - 1\right)
 *  \right],
 * @f]
 * with @f$p_1,p_2@f$ free parameters and @f$C@f$ a normalisation constant:
 * @f[
 *  C = \left( \int_0^\pi{} \left[
 *    p_1 + \frac{5}{4} p_2 \left(\cos^2(\beta{}) - 1\right)
 *  \right] \sin(\beta{}) d\beta{} \right)^{-1} =
 *  \left( 2 p_1 - \frac{5}{3} p_2 \right)^{-1}.
 * @f]
 *
 * The free parameter @f$p_2@f$ is linked to the distribution parameter
 * @f[
 *  \langle{} \cos^2(\beta{}) \rangle{} =
 *    \int_0^\pi{} \sin(\beta{}) \cos^2(\beta{}) p(\beta{}) d\beta{} =
 *    \frac{1}{3}C (2p_1-p_2)
 * @f]
 * (often quoted in literature). The parameter @f$p_1@f$ is determined from
 * the requirement
 * @f[
 *  p(0) = p(\pi{}) = Cp_1 = 0
 * @f]
 * for prolate spheroids, and
 * @f[
 *  p\left(\frac{\pi{}}{2}\right) = C\left(p_1 - \frac{5}{4}p_2\right) = 0
 * @f]
 * for oblate spheroids.
 *
 * As in Mishchenko (1991), we will parameterise the distribution function
 * using @f$\langle{} \cos^2(\beta{}) \rangle{}@f$.
 *
 * Note that Mishchenko (1991) quotes the wrong expression for
 * @f$\langle{} \cos^2(\beta{}) \rangle{}@f$ and sets the parameter
 * @f$p_1=\frac{1}{2}@f$ for both oblate and prolate spheroids. The results in
 * Table 1 in this paper can however only be reproduced using the general
 * expression given above.
 *
 * Since OrientationDistribution computes the normalisation constant @f$C@f$
 * internally, we can set it to an arbitrary value and only need to compute the
 * parameters @f$p_1@f$ and @f$p_2@f$. For oblate spheroids, these are given by
 * @f[
 *   p_1 = \frac{5}{2} \langle{} \cos^2(\beta{}) \rangle{},
 * @f]
 * @f[
 *   p_2 = 2 \langle{} \cos^2(\beta{}) \rangle{}.
 * @f]
 * For prolate spheroids, we have
 * @f[
 *   p_1 = 0,
 * @f]
 * @f[
 *   p_2 = -3 \langle{} \cos^2(\beta{}) \rangle{}.
 * @f]
 */
class MishchenkoOrientationDistribution : public OrientationDistribution {
private:
  /*! @brief Parameter @f$p_1@f$. */
  const float_type _p1;

  /*! @brief Parameter @f$p_2@f$. */
  const float_type _p2;

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

    return _p1 - 1.25 * _p2 * sinbeta * sinbeta;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   * @param axis_ratio Axis ratio @f$d@f$ of the spheroidal scatterers, used to
   * determine whether the prolate (@f$d < 1@f$) or oblate (@f$d > 1@f$)
   * parameter formulae should be used.
   * @param cos2beta Parameter @f$\langle{} \cos^2(\beta{}) \rangle{}@f$.
   */
  inline MishchenkoOrientationDistribution(const uint_fast32_t nmax,
                                           const float_type axis_ratio,
                                           const float_type cos2beta)
      : OrientationDistribution(nmax),
        _p1((axis_ratio < 1.) ? 0. : 2.5 * cos2beta),
        _p2((axis_ratio < 1.) ? -3. * cos2beta : 2. * cos2beta) {

    initialise();
  }

  virtual ~MishchenkoOrientationDistribution() {}
};

#endif // MISHCHENKOORIENTATIONDISTRIBUTION_HPP
