/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Matrix.hpp
 *
 * @brief Matrix functionality.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Error.hpp"

#include <cinttypes>
#include <complex>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Matrix class.
 *
 * Used to represent a 2D matrix of fixed dimensions (dimensions can be set
 * during construction).
 *
 * @tparam DATA_TYPE Data type of the matrix elements.
 */
template <typename DATA_TYPE> class Matrix {
private:
  /*! @brief Number of rows in the matrix. */
  const uint_fast32_t _number_of_rows;

  /*! @brief Number of columns in the matrix. */
  const uint_fast32_t _number_of_columns;

  /*! @brief Flat array containing the matrix elements. */
  std::vector<DATA_TYPE> _array;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_rows Number of rows in the matrix.
   * @param number_of_columns Number of columns in the matrix.
   */
  inline Matrix(const uint_fast32_t number_of_rows,
                const uint_fast32_t number_of_columns)
      : _number_of_rows(number_of_rows), _number_of_columns(number_of_columns) {
    _array.resize(_number_of_rows * _number_of_columns, DATA_TYPE(0.));
  }

  /**
   * @brief Clear the contents of the matrix.
   */
  inline void reset() {
    for (uint_fast32_t i = 0; i < _array.size(); ++i) {
      _array[i] = 0.;
    }
  }

  /**
   * @brief Get the number of rows in the matrix.
   *
   * @return Number of rows in the matrix.
   */
  inline uint_fast32_t get_number_of_rows() const { return _number_of_rows; }

  /**
   * @brief Get the number of columns in the matrix.
   *
   * @return Number of columns in the matrix.
   */
  inline uint_fast32_t get_number_of_columns() const {
    return _number_of_columns;
  }

  /**
   * @brief Access operator.
   *
   * This function performs the bookkeeping required to figure out where in the
   * flat data array the element is stored.
   *
   * We store the elements per row, i.e. first all columns of the first row,
   * then all the columns of the second row...
   *
   * Note that if all access to the matrix elements goes through this function,
   * we can in principle change this order very quickly.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the corresponding matrix element.
   */
  inline DATA_TYPE &operator()(const uint_fast32_t i, const uint_fast32_t j) {
    return _array[i * _number_of_columns + j];
  }

  /**
   * @brief Access operator.
   *
   * This function performs the bookkeeping required to figure out where in the
   * flat data array the element is stored.
   *
   * We store the elements per row, i.e. first all columns of the first row,
   * then all the columns of the second row...
   *
   * Note that if all access to the matrix elements goes through this function,
   * we can in principle change this order very quickly.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the corresponding matrix element (read only).
   */
  inline const DATA_TYPE &operator()(const uint_fast32_t i,
                                     const uint_fast32_t j) const {
    return _array[i * _number_of_columns + j];
  }

  /**
   * @brief Get the row with the given index as a pointer that can be accessed
   * like an array.
   *
   * @param i Row index.
   * @return Pointer to the beginning of that row.
   */
  inline DATA_TYPE *get_row(const uint_fast32_t i) {
    return &_array[i * _number_of_columns];
  }

  /**
   * @brief Dump the Matrix to a binary file with the given name.
   *
   * The first two 32-bit values in the output file will contain the number of
   * rows and columns in the matrix, followed by the size of a single element
   * of the matrix, followed by all the matrix elements in C++ array order (i.e.
   * first all columns of the first row, and so on...).
   *
   * @param filename Name of the output file.
   */
  inline void binary_dump(const std::string filename) const {

    std::ofstream file(filename, std::ofstream::binary);
    // make sure these numbers are actual 32-bit integers before writing them
    // out
    const uint32_t number_of_rows = _number_of_rows;
    file.write(reinterpret_cast<const char *>(&number_of_rows),
               sizeof(number_of_rows));
    const uint32_t number_of_columns = _number_of_columns;
    file.write(reinterpret_cast<const char *>(&number_of_columns),
               sizeof(number_of_columns));
    const uint32_t datasize = sizeof(DATA_TYPE);
    file.write(reinterpret_cast<const char *>(&datasize), sizeof(datasize));
    const uint_fast32_t size = _number_of_columns * _number_of_rows;
    file.write(reinterpret_cast<const char *>(&_array[0]), size * datasize);
  }

  /**
   * @brief Invert the matrix using an auxiliary PLU decomposition.
   *
   * This only works if the matrix is square. If not, an error is thrown.
   *
   * This is a 4 step algorithm. For what follows, we will use the symbol
   * @f$A@f$ to represent the matrix, and @f$A_{ij}@f$ (@f$i,j \in{} [0, N[@f$),
   * with @f$N@f$ the size of the matrix in one of its dimensions).
   *
   * First we decompose the matrix into three matrices:
   * @f[
   *    A = P \times{} L \times{} U,
   * @f]
   * where @f$P@f$ is a permutation matrix that simply interchanges two rows
   * when it is premultiplied with another matrix, @f$L@f$ is a lower
   * triangular matrix with ones along its diagonal, and @f$U@f$ is an upper
   * triangular matrix with no restrictions on its diagonal elements. Any
   * regular square matrix can be decomposed in this form. This is called
   * an LU decomposition with partial pivoting.
   *
   * To perform the decomposition, we use an in-place version of the Doolittle
   * algorithm
   * (https://en.wikipedia.org/wiki/LU_decomposition#Doolittle_algorithm),
   * where we perform a single loop over all rows, and for every @f$i@f$th row
   * row select the @f$j@f$th row (@f$j \in{} [i, N[@f$) that has the element in
   * the @f$i@f$th column with the highest absolute value (this is the pivot).
   * We then interchange rows @f$i@f$ and @f$j@f$, and store the index @f$j@f$
   * in the permutation matrix (represented as a pivot array of size @f$N@f$).
   * We then divide all elements in the @f$i@f$th column in rows @f$j>i@f$
   * by the element @f$A_{ii}@f$ (these are the values in the matrix @f$L@f$
   * for this column), and eliminate the elements below the diagonal from
   * the original matrix by subtracting @f$A_{ji}/A_{ii}@f$ times the original
   * row from all rows @f$j>i@f$ (we only do this for the columns @f$k>i@f$ so
   * that the elements of the @f$L@f$ matrix can be stored in the lower part
   * of the original matrix.
   *
   * At the end of step 1, we are left with a modified matrix @f$A@f$ that now
   * contains the elements of @f$L@f$ in its lower part, and the elements of
   * @f$U@f$ in the upper part (including the diagonal).
   *
   * To compute the inverse of the original matrix, we make use of
   * @f[
   *    A^{-1} = (P \times{} L \times{} u)^{-1},
   * @f]
   * and rewrite this as
   * @f[
   *    A^{-1} \times{} P \times{} L = (P^{-1} \times{} A)^{-1} \times{} L =
   *      U^{-1}.
   * @f]
   *
   * In step 2, we first invert the upper triangular matrix @f$U@f$, using
   * a forward substitution algorithm. The diagonal elements of @f$U^{-1}@f$
   * are given by
   * @f[
   *    U^{-1}_{ii} = \frac{1}{U_{ii}},
   * @f]
   * while the upper triangular elements are
   * @f[
   *    U^{-1}_{ij} = -\frac{1}{U_{jj}} \sum_{k=i}^j U^{-1}_{ik} U_{kj}.
   * @f]
   *
   * We then solve the equation
   * @f[
   *    (P^{-1} \times{} A)^{-1} \times{} L = U^{-1},
   * @f]
   * using a backward substitution algorithm:
   * @f[
   *    (P^{-1} \times{} A)^{-1}_{ij} = U^{-1}_{ij} - \sum_{k=j+1}^{N}
   *      (P^{-1} \times{} A)^{-1}_{ik} L_{kj}.
   * @f]
   * This is step 3.
   *
   * In step 4, we solve
   * @f[
   *    A^{-1} = (P^{-1} \times{} A)^{-1} \times{} P^{-1}.
   * @f]
   * Due to the simplicity of the permutation matrix (it is a permutation of
   * the unit matrix), this boils down to simply swapping the columns of the
   * matrix in reverse order of the pivot array.
   *
   * @param number_of_columns Number of columns (and rows) in the current
   * matrix. This does not need to match the actual number of rows and columns
   * in the data array, to allow for the inversion of smaller matrices.
   * @param pivot_array Preallocated array to store pivot values in.
   * @param pivot_array_size Size of the pivot array.
   * @param work Preallocated array to store temporary columns in.
   * @param work_size Size of the work array.
   */
  inline void plu_inverse(const uint_fast32_t number_of_columns,
                          uint_fast32_t *pivot_array,
                          const uint_fast32_t pivot_array_size, DATA_TYPE *work,
                          const uint_fast32_t work_size) {

    ctm_assert(number_of_columns <= _number_of_columns);
    ctm_assert(number_of_columns <= _number_of_rows);
    ctm_assert(pivot_array_size >= number_of_columns);
    ctm_assert(work_size >= number_of_columns);

    // required so that the compiler finds the correct division operator for
    // std::complex<boost::multiprecision::cpp_bin_float_quad> DATA_TYPEs
    const DATA_TYPE one(1.);

    // alias the object using the label A (to be consistent with the notation
    // in the function documentation)
    Matrix &A = *this;

    // step 1: PLU decomposition of the original matrix
    for (uint_fast32_t i = 0; i < number_of_columns; ++i) {
      // find the next pivot
      uint_fast32_t imax = i;
      // declaring the variables Smax and this_Smax as auto is the only way of
      // making sure that we use the right type for
      // std::complex<boost::multiprecision::cpp_bin_float_quad>> DATA_TYPEs
      auto Smax = abs(A(imax, i));
      for (uint_fast32_t j = i + 1; j < number_of_columns; ++j) {
        const auto this_Smax = abs(A(j, i));
        if (this_Smax > Smax) {
          Smax = this_Smax;
          imax = j;
        }
      }
      // store the pivot in the pivot array
      pivot_array[i] = imax;
      // swap rows if necessary
      if (i != imax) {
        for (uint_fast32_t j = 0; j < number_of_columns; ++j) {
          // we need a temporary variable to now overwrite the original A_ij
          const DATA_TYPE temp = A(i, j);
          A(i, j) = A(imax, j);
          A(imax, j) = temp;
        }
      }
      // compute the inverse of A_ii to save on divisions
      const DATA_TYPE Aii_inv = one / A(i, i);
      // now multiply all elements below row i in column i with this value
      for (uint_fast32_t j = i + 1; j < number_of_columns; ++j) {
        A(j, i) *= Aii_inv;
      }
      // eliminate lower triangular elements in column i from the future U
      // matrix
      for (uint_fast32_t j = i + 1; j < number_of_columns; ++j) {
        // since we store the L matrix in the lower diagonal part of the matrix,
        // we only update columns k > i
        for (uint_fast32_t k = i + 1; k < number_of_columns; ++k) {
          A(j, k) -= A(j, i) * A(i, k);
        }
      }
    }

    // step 2: inversion of the U part of A
    for (uint_fast32_t i = 0; i < number_of_columns; ++i) {
      // compute the diagonal element of the inverse matrix for row i
      A(i, i) = one / A(i, i);
      // now use forward substitution to compute the other elements in this row
      for (uint_fast32_t j = i + 1; j < number_of_columns; ++j) {
        // we need a temporary variable, as A_ij also features in the loop below
        DATA_TYPE Aij(0.);
        for (uint_fast32_t k = i; k < j; ++k) {
          Aij -= A(i, k) * A(k, j);
        }
        // we cannot save on divisions here, as we still need the original
        // A_jj later on
        A(i, j) = Aij / A(j, j);
      }
    }

    // step 3: solve inv(A)*L = inv(U)
    // we need a temporary array to store the new columns while we update them
    // note that we need to iterate backwards for this algorithm
    for (uint_fast32_t jp1 = number_of_columns; jp1 > 0; --jp1) {
      // unsigned counters overflow when you decrement them below 0, so we
      // need a second variable to handle j = 0 in the loop condition
      const uint_fast32_t j = jp1 - 1;
      // initialise the new column with the values in the old column above the
      // diagonal
      for (uint_fast32_t i = 0; i < jp1; ++i) {
        work[i] = A(i, j);
      }
      // zero out the other elements in the new column
      for (uint_fast32_t i = jp1; i < number_of_columns; ++i) {
        work[i] = 0.;
      }
      // now apply the summation to eliminate unknowns in the column
      for (uint_fast32_t k = jp1; k < number_of_columns; ++k) {
        for (uint_fast32_t i = 0; i < number_of_columns; ++i) {
          work[i] -= A(i, k) * A(k, j);
        }
      }
      // we are done with the old column, overwrite it with the new column
      for (uint_fast32_t i = 0; i < number_of_columns; ++i) {
        A(i, j) = work[i];
      }
    }

    // step 4: exchange columns based on the pivot array
    // again, this is a backward iteration
    // note that we swap columns, since we multiply the inverse permutation
    // matrix on the other side
    for (uint_fast32_t jp1 = number_of_columns; jp1 > 0; --jp1) {
      // again, we need a second variable to handle unsigned integer overflow
      const uint_fast32_t j = jp1 - 1;
      const uint_fast32_t jp = pivot_array[j];
      // only swap columns if the pivot index is different
      if (jp != j) {
        for (uint_fast32_t i = 0; i < number_of_columns; ++i) {
          const DATA_TYPE temp = A(i, j);
          A(i, j) = A(i, jp);
          A(i, jp) = temp;
        }
      }
    }
  }

  /**
   * @brief PLU inversion function that uses its own pivot and work arrays.
   *
   * @param number_of_columns Number of columns (and rows) in the current
   * matrix. This does not need to match the actual number of rows and columns
   * in the data array, to allow for the inversion of smaller matrices.
   */
  inline void plu_inverse(const uint_fast32_t number_of_columns) {

    std::vector<uint_fast32_t> pivot_array(number_of_columns);
    std::vector<DATA_TYPE> work(number_of_columns);
    plu_inverse(number_of_columns, &pivot_array[0], number_of_columns, &work[0],
                number_of_columns);
  }
};

#endif // MATRIX_HPP
