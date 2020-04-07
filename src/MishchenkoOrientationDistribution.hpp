/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

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
 *  p(\beta{}) =
 *    \frac{1}{2} + \frac{5}{4} p_2 \left(3\cos^2(\beta{}) - 1\right),
 * @f]
 * with @f$p_2@f$ a free parameter. This distribution satisfies the
 * normalisation condition
 * @f[
 *  \int_0^{\pi{}} p(\beta{}) \sin(\beta{}) d\beta{} = 1,
 * @f]
 * as can be easily checked (note that the expression given in Mishchenko (1991)
 * is missing the factor @f$3@f$, while quoting the right expression in terms
 * of the Legendre polynomial @f$P_2(\cos(\beta{}))@f$; this is likely a typo).
 *
 * The free parameter @f$p_2@f$ is linked to the distribution parameter
 * @f[
 *  \langle{} \cos^2(\beta{}) \rangle{} =
 *    \int_0^\pi{} \sin(\beta{}) \cos^2(\beta{}) p(\beta{}) d\beta{} =
 *    \frac{1}{3} (2p_2 + 1)
 * @f]
 * (often quoted in literature).
 *
 * As in Mishchenko (1991), we will parameterise the distribution function
 * using @f$\langle{} \cos^2(\beta{}) \rangle{}@f$ and compute the distribution
 * parameter from
 * @f[
 *  p_2 = \frac{3}{2} \langle{} \cos^2(\beta{}) \rangle{} - \frac{1}{2}.
 * @f]
 */
class MishchenkoOrientationDistribution : public OrientationDistribution {
private:
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

    return 0.5 + 1.25 * _p2 * (3. * cosbeta * cosbeta - 1.);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   * @param cos2beta Parameter @f$\langle{} \cos^2(\beta{}) \rangle{}@f$.
   */
  inline MishchenkoOrientationDistribution(const uint_fast32_t nmax,
                                           const float_type cos2beta)
      : OrientationDistribution(nmax), _p2(1.5 * cos2beta - 0.5) {

    ctm_assert(cos2beta >= 0.2);
    ctm_assert(cos2beta <= 0.6);

    // note that we choose not to compute the coefficients using the integral,
    // since we already know what the values of the coefficients should be...
    for (uint_fast32_t i = 0; i < _coefficients.size(); ++i) {
      _coefficients[i] = 0.;
    }
    _coefficients[0] = 1.;
    _coefficients[2] = _p2;
  }

  virtual ~MishchenkoOrientationDistribution() {}
};

#endif // MISHCHENKOORIENTATIONDISTRIBUTION_HPP
