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
  /*! @brief Silicon properties table. */
  LinearInterpolatedTable<5, float_type> _table_silicon;

  /*! @brief Table for small, parallel aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_parallel_small;

  /*! @brief Table for large, parallel aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_parallel_large;

  /*! @brief Table for small, perpendicularly aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_perpendicular_small;

  /*! @brief Table for large, perpendicularly aligned carbon grains. */
  LinearInterpolatedTable<5, float_type> _table_carbon_perpendicular_large;

public:
  /**
   * @brief Constructor.
   */
  inline DraineDustProperties() {

    _table_silicon.from_ascii_file(DRAINEDUSTPROPERTIESDATALOCATION
                                   "callindex.out_silD03");
    // convert wavelength from micron to m
    _table_silicon.multiply_column<0>(1.e-6);
    _table_silicon.add_column<1>(1.);
    _table_silicon.add_column<3>(1.);

    _table_carbon_parallel_small.from_ascii_file(
        DRAINEDUSTPROPERTIESDATALOCATION "callindex.out_CpaD03_0.01");
    // convert wavelength from micron to m
    _table_carbon_parallel_small.multiply_column<0>(1.e-6);
    _table_carbon_parallel_small.add_column<1>(1.);
    _table_carbon_parallel_small.add_column<3>(1.);

    _table_carbon_parallel_large.from_ascii_file(
        DRAINEDUSTPROPERTIESDATALOCATION "callindex.out_CpaD03_0.10");
    // convert wavelength from micron to m
    _table_carbon_parallel_large.multiply_column<0>(1.e-6);
    _table_carbon_parallel_large.add_column<1>(1.);
    _table_carbon_parallel_large.add_column<3>(1.);

    _table_carbon_perpendicular_small.from_ascii_file(
        DRAINEDUSTPROPERTIESDATALOCATION "callindex.out_CpeD03_0.01");
    // convert wavelength from micron to m
    _table_carbon_perpendicular_small.multiply_column<0>(1.e-6);
    _table_carbon_perpendicular_small.add_column<1>(1.);
    _table_carbon_perpendicular_small.add_column<3>(1.);

    _table_carbon_perpendicular_large.from_ascii_file(
        DRAINEDUSTPROPERTIESDATALOCATION "callindex.out_CpeD03_0.10");
    // convert wavelength from micron to m
    _table_carbon_perpendicular_large.multiply_column<0>(1.e-6);
    _table_carbon_perpendicular_large.add_column<1>(1.);
    _table_carbon_perpendicular_large.add_column<3>(1.);
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
      std::complex<float_type> small_value, large_value;
      if (grain_size < 1.e-7) {
        const TableRow<5, float_type> &row =
            _table_carbon_parallel_small.get_row<0>(wavelength);
        small_value = std::complex<float_type>(row[3], row[4]);
      }
      if (grain_size > 1.e-8) {
        const TableRow<5, float_type> &row =
            _table_carbon_parallel_large.get_row<0>(wavelength);
        large_value = std::complex<float_type>(row[3], row[4]);
      }
      float_type small_fraction, large_fraction;
      if (grain_size <= 1.e-8) {
        small_fraction = 1.;
        large_fraction = 0.;
      } else if (grain_size < 1.e-7) {
        small_fraction = (grain_size * 1.e6 - 0.01) / 0.09;
        large_fraction = 1. - small_fraction;
      } else {
        small_fraction = 0.;
        large_fraction = 1.;
      }
      return small_fraction * small_value + large_fraction + large_value;
    } else if (grain_type == DUSTGRAINTYPE_CARBON_PERPENDICULAR) {
      std::complex<float_type> small_value, large_value;
      if (grain_size < 1.e-7) {
        const TableRow<5, float_type> &row =
            _table_carbon_perpendicular_small.get_row<0>(wavelength);
        small_value = std::complex<float_type>(row[3], row[4]);
      }
      if (grain_size > 1.e-8) {
        const TableRow<5, float_type> &row =
            _table_carbon_perpendicular_large.get_row<0>(wavelength);
        large_value = std::complex<float_type>(row[3], row[4]);
      }
      float_type small_fraction, large_fraction;
      if (grain_size <= 1.e-8) {
        small_fraction = 1.;
        large_fraction = 0.;
      } else if (grain_size < 1.e-7) {
        small_fraction = (grain_size * 1.e6 - 0.01) / 0.09;
        large_fraction = 1. - small_fraction;
      } else {
        small_fraction = 0.;
        large_fraction = 1.;
      }
      return small_fraction * small_value + large_fraction + large_value;
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
