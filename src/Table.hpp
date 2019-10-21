/**
 * @file Table.hpp
 *
 * @brief Table.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TABLE_HPP
#define TABLE_HPP

#include "Error.hpp"

#include <cinttypes>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Row of the table.
 *
 * @tparam NUMBER_OF_COLUMNS Number of columns in each row.
 * @tparam DATA_TYPE Data type of the row elements.
 */
template <uint_fast8_t NUMBER_OF_COLUMNS, typename DATA_TYPE> class TableRow {
private:
  /*! @brief Columns of the row. */
  DATA_TYPE _data[NUMBER_OF_COLUMNS];

public:
  /**
   * @brief Empty constructor.
   */
  inline TableRow() {
    for (uint_fast8_t i = 0; i < NUMBER_OF_COLUMNS; ++i) {
      _data[i] = 0.;
    }
  }

  /**
   * @brief Access an element in the row.
   *
   * @param column Column index.
   * @return Reference to the corresponding row element.
   */
  inline DATA_TYPE &operator[](const uint_fast8_t column) {
    ctm_assert(column < NUMBER_OF_COLUMNS);
    return _data[column];
  }

  /**
   * @brief Access an element in the row.
   *
   * @param column Column index.
   * @return Reference to the corresponding row element.
   */
  inline const DATA_TYPE &operator[](const uint_fast8_t column) const {
    ctm_assert(column < NUMBER_OF_COLUMNS);
    return _data[column];
  }
};

/**
 * @brief Table class.
 *
 * Used to represent a table with a template number of columns, adjustable
 * number of rows and template cell data type.
 *
 * @tparam NUMBER_OF_COLUMNS Number of columns in each row.
 * @tparam DATA_TYPE Data type of the table elements.
 */
template <uint_fast8_t NUMBER_OF_COLUMNS, typename DATA_TYPE> class Table {
protected:
  /*! @brief Table rows. */
  std::vector<TableRow<NUMBER_OF_COLUMNS, DATA_TYPE>> _rows;

public:
  /**
   * @brief Access the row with the given index.
   *
   * @param row Row index.
   * @return Corresponding row of the table.
   */
  inline TableRow<NUMBER_OF_COLUMNS, DATA_TYPE> &operator[](const size_t row) {
    ctm_assert(row < _rows.size());
    return _rows[row];
  }

  /**
   * @brief Access the row with the given index.
   *
   * @param row Row index.
   * @return Corresponding row of the table.
   */
  inline const TableRow<NUMBER_OF_COLUMNS, DATA_TYPE> &
  operator[](const size_t row) const {
    ctm_assert(row < _rows.size());
    return _rows[row];
  }

  /**
   * @brief Fill the Table with the contents of the ASCII file with the given
   * name.
   *
   * @param filename Name of the file.
   */
  inline void from_ascii_file(const std::string filename) {
    std::ifstream ifile(filename);
    if (!ifile.good()) {
      ctm_error("Unable to open file: \"%s\"!", filename.c_str());
    }
    std::string line;
    while (std::getline(ifile, line)) {
      // skip comment lines
      if (line[0] != '#') {
        std::istringstream linestr(line);
        _rows.resize(_rows.size() + 1);
        for (uint_fast8_t i = 0; i < NUMBER_OF_COLUMNS; ++i) {
          linestr >> _rows.back()[i];
        }
      }
    }
  }

  /**
   * @brief Get the number of rows in the table.
   *
   * @return Number of rows in the table.
   */
  inline virtual size_t size() const { return _rows.size(); }

  /**
   * @brief Multiply all elements in a column with the given factor.
   *
   * @param factor Multiplicative factor.
   * @tparam COLUMN_NUMBER Column that needs to be multiplied.
   */
  template <uint_fast8_t COLUMN_NUMBER>
  inline void multiply_column(const DATA_TYPE factor) {
    for (uint_fast32_t i = 0; i < _rows.size(); ++i) {
      _rows[i][COLUMN_NUMBER] *= factor;
    }
  }
};

/**
 * @brief Linearly interpolated table class.
 *
 * Used to represent a table with a template number of columns, adjustable
 * number of rows and template cell data type, for which elements can be
 * accessed using linear interpolation.
 *
 * @tparam NUMBER_OF_COLUMNS Number of columns in each row.
 * @tparam DATA_TYPE Data type of the table elements.
 */
template <uint_fast8_t NUMBER_OF_COLUMNS, typename DATA_TYPE>
class LinearInterpolatedTable : public Table<NUMBER_OF_COLUMNS, DATA_TYPE> {
public:
  /**
   * @brief Get the row corresponding to the given x value, using linear
   * interpolation to determine the y values.
   *
   * @param xvalue X value.
   * @return Row containing the corresponding y values.
   * @tparam INTERPOLATING_COLUMN Index of the column that is used for the x
   * values.
   */
  template <uint_fast8_t INTERPOLATING_COLUMN>
  inline TableRow<NUMBER_OF_COLUMNS, DATA_TYPE>
  get_row(const DATA_TYPE xvalue) const {

    size_t ilow = 0;
    size_t ihigh = Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows.size() - 1;
    DATA_TYPE xlow =
        Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[ilow][INTERPOLATING_COLUMN];
    DATA_TYPE xhigh =
        Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[ihigh][INTERPOLATING_COLUMN];
    if (xlow < xhigh) {
      while (ilow + 1 != ihigh) {
        const size_t inext = (ilow + ihigh) >> 1;
        const DATA_TYPE xnext =
            Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[inext]
                                                      [INTERPOLATING_COLUMN];
        if (xnext <= xvalue) {
          xlow = xnext;
          ilow = inext;
        } else {
          xhigh = xnext;
          ihigh = inext;
        }
      }
    } else if (xlow > xhigh) {
      while (ilow + 1 != ihigh) {
        const size_t inext = (ilow + ihigh) >> 1;
        const DATA_TYPE xnext =
            Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[inext]
                                                      [INTERPOLATING_COLUMN];
        if (xnext >= xvalue) {
          xlow = xnext;
          ilow = inext;
        } else {
          xhigh = xnext;
          ihigh = inext;
        }
      }
    } else {
      ctm_error("Table endpoints are the same!");
    }

    const float_type xfac = (xvalue - xlow) / (xhigh - xlow);
    TableRow<NUMBER_OF_COLUMNS, DATA_TYPE> result;
    for (uint_fast8_t i = 0; i < NUMBER_OF_COLUMNS; ++i) {
      result[i] =
          xfac * Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[ilow][i] +
          (1. - xfac) * Table<NUMBER_OF_COLUMNS, DATA_TYPE>::_rows[ihigh][i];
    }
    return result;
  }
};

#endif // TABLE_HPP
