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
 * @file DraineDustProperties.hpp
 *
 * @brief Draine and collaborators dust optical dust properties.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef DRAINEDUSTPROPERTIES_HPP
#define DRAINEDUSTPROPERTIES_HPP

#include "DraineDustPropertiesDataLocation.hpp"
#include "DustProperties.hpp"
#include "Table.hpp"

/**
 * @brief Draine and collaborators dust optical dust properties.
 */
class DraineDustProperties : public DustProperties {
private:
  /*! @brief Temperature of the dust (in K). */
  const float_type _dust_temperature;

  /*! @brief Silicon properties table. */
  LinearInterpolatedTable<5, float_type> _table_silicon;

  /*! @brief Table for parallel aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_parallel;

  /*! @brief Table for perpendicularly aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_perpendicular;

  /**
   * @brief Get the dielectric function contribution from free electrons to
   * the carbon dielectric function.
   *
   * The correction is based on equation 2.2 in Draine & Lee (1984), using the
   * parameter values tabulated in their table 1 (and noting the typo mentioned
   * in the erratum for that paper).
   *
   * The grain size dependence is as given in equations 2.5-2.6. The resonances
   * in equation 2.7 are independent of temperature and grain size and are hence
   * not part our correction.
   *
   * @param T Temperature (in K).
   * @param wavelength Wavelength at which the dielectric function is evaluated
   * (in m).
   * @param grain_size Size of the dust grain (in m).
   * @param is_parallel Is the dielectric function evaluated for radiation that
   * falls in parallel to the "c-axis"?
   * @return Dielectric function contribution from free electrons.
   */
  inline static std::complex<float_type>
  get_C_epsilon_free(const float_type T, const float_type wavelength,
                     const float_type grain_size, const bool is_parallel) {

    // Draine & Lee (1984), table 1
    float_type omega_p, tau_bulk, TF, vF;
    if (is_parallel) {
      omega_p = 1.53e14;
      tau_bulk = 1.4e-14;
      TF = 255.;
      // we converted from cm s^-1 to m s^-1
      vF = 3.7e4;
    } else {
      const float_type T2 = T * T;
      // Erratum, 1987, ApJ, 318, 485
      omega_p = 4.33e14 * sqrt(1. - 6.24e-3 * T + 3.66e-5 * T2);
      tau_bulk = 4.2e-11 / (1. + 0.322 * T + 0.0013 * T2);
      TF = 255.;
      vF = 4.5e5;
    }

    // size-based correction for tau
    const float_type veff = vF * sqrt(1. + T / TF);
    tau_bulk = 1. / (1. / tau_bulk + veff / grain_size);

    // c = 299792458 m s^-1
    const float_type omega = 2. * M_PI * 299792458. / wavelength;
    const float_type omega_tau = omega * tau_bulk;
    const float_type omega_p_tau = omega_p * tau_bulk;
    // constructing complex values from doubles is not straightforward in C++,
    // so this is the cleanest way I could find.
    const std::complex<float_type> epsilon_corr_nom(-omega_p_tau * omega_p_tau,
                                                    0.);
    const std::complex<float_type> epsilon_corr_denom(omega_tau * omega_tau,
                                                      omega_tau);
    return epsilon_corr_nom / epsilon_corr_denom;
  }

  /**
   * @brief Convert the given dielectric function to the equivalent refractive
   * index.
   *
   * If the dielectric function is given by
   * @f[
   *  \varepsilon{} = \varepsilon{}_r + i \varepsilon{}_i,
   * @f]
   * then the refractive index is
   * @f[
   *  n = n_r + i n_i,
   * @f]
   * with
   * @f[
   *  n_r = \sqrt{\frac{|\varepsilon{}| + \varepsilon{}_r}{2}},
   * @f]
   * @f[
   *  n_i = \sqrt{\frac{|\varepsilon{}| - \varepsilon{}_r}{2}},
   * @f]
   * and
   * @f[
   *  |\varepsilon{}| = \sqrt{\varepsilon{}_r^2 + \varepsilon{}_i^2}.
   * @f]
   *
   * @param eps Input dielectric function.
   * @return Equivalent refractive index.
   */
  inline static std::complex<float_type>
  eps_to_mr(const std::complex<float_type> eps) {
    const float_type epsr = eps.real();
    const float_type epsi = eps.imag();
    const float_type epsabs = sqrt(epsr * epsr + epsi * epsi);
    return std::complex<float_type>(sqrt(0.5 * (epsabs + epsr)),
                                    sqrt(0.5 * (epsabs - epsr)));
  }

public:
  /**
   * @brief Constructor.
   *
   * @param dust_temperature Temperature of the dust (in K).
   */
  inline DraineDustProperties(const float_type dust_temperature = 20.)
      : _dust_temperature(dust_temperature) {

    _table_silicon.from_ascii_file(DRAINEDUSTPROPERTIESDATALOCATION
                                   "callindex.out_silD03");
    // convert wavelength from micron to m
    _table_silicon.multiply_column<0>(1.e-6);
    _table_silicon.add_column<1>(1.);
    _table_silicon.add_column<3>(1.);

    _table_carbon_parallel.from_ascii_file(DRAINEDUSTPROPERTIESDATALOCATION
                                           "callindex.out_CpaD03_0.10");
    // convert wavelength from micron to m
    _table_carbon_parallel.multiply_column<0>(1.e-6);
    _table_carbon_parallel.add_column<1>(1.);
    _table_carbon_parallel.add_column<3>(1.);

    // remove the temperature and size dependent contribution to the dielectric
    // function
    // the values have been tabulated for 20 K and a grain size of 0.1 micron
    // and include the correction term for those parameter values
    for (size_t i = 0; i < _table_carbon_parallel.size(); ++i) {
      const float_type wavelength = _table_carbon_parallel[i][0];
      const std::complex<float_type> epsilon_corr =
          get_C_epsilon_free(20., wavelength, 1.e-7, true);
      _table_carbon_parallel[i][1] -= epsilon_corr.real();
      _table_carbon_parallel[i][2] -= epsilon_corr.imag();
    }

    _table_carbon_perpendicular.from_ascii_file(DRAINEDUSTPROPERTIESDATALOCATION
                                                "callindex.out_CpeD03_0.10");
    // convert wavelength from micron to m
    _table_carbon_perpendicular.multiply_column<0>(1.e-6);
    _table_carbon_perpendicular.add_column<1>(1.);
    _table_carbon_perpendicular.add_column<3>(1.);

    // remove the temperature and size dependent contribution to the dielectric
    // function
    // the values have been tabulated for 20 K and a grain size of 0.1 micron
    // and include the correction term for those parameter values
    for (size_t i = 0; i < _table_carbon_perpendicular.size(); ++i) {
      const float_type wavelength = _table_carbon_perpendicular[i][0];
      const std::complex<float_type> epsilon_corr =
          get_C_epsilon_free(20., wavelength, 1.e-7, false);
      _table_carbon_perpendicular[i][1] -= epsilon_corr.real();
      _table_carbon_perpendicular[i][2] -= epsilon_corr.imag();
    }
  }

  virtual ~DraineDustProperties() {}

  /**
   * @brief Get the refractive index for the given wavelength, grain size and
   * grain type.
   *
   * @param wavelength Wavelength of the incoming radiation (in m).
   * @param grain_size Size of the grain (in m).
   * @param grain_type Type of dust grain.
   * @return Complex refractive index of the grain for radiation at this
   * wavelength.
   */
  inline std::complex<float_type>
  get_refractive_index(const float_type wavelength, const float_type grain_size,
                       const int_fast32_t grain_type) const {

    if (grain_type == DUSTGRAINTYPE_CARBON_PARALLEL) {
      const TableRow<5, float_type> &row =
          _table_carbon_parallel.get_row<0>(wavelength);
      std::complex<float_type> eps(row[1], row[2]);
      const std::complex<float_type> epscorr =
          get_C_epsilon_free(_dust_temperature, wavelength, grain_size, true);
      eps += epscorr;
      return eps_to_mr(eps);
    } else if (grain_type == DUSTGRAINTYPE_CARBON_PERPENDICULAR) {
      const TableRow<5, float_type> &row =
          _table_carbon_perpendicular.get_row<0>(wavelength);
      std::complex<float_type> eps(row[1], row[2]);
      const std::complex<float_type> epscorr =
          get_C_epsilon_free(_dust_temperature, wavelength, grain_size, false);
      eps += epscorr;
      return eps_to_mr(eps);
    } else if (grain_type == DUSTGRAINTYPE_SILICON) {
      const TableRow<5, float_type> &row =
          _table_silicon.get_row<0>(wavelength);
      return std::complex<float_type>(row[3], row[4]);
    } else {
      ctm_error("Unknown dust grain type: %" PRIiFAST32 "!", grain_type);
      return std::complex<float_type>(0.);
    }
  }
};

#endif // DRAINEDUSTPROPERTIES_HPP
