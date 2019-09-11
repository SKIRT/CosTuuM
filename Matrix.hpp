/**
 * @file Matrix.hpp
 *
 * @brief Matrix functionality.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include <cinttypes>
#include <complex>
#include <vector>

/**
 * @brief Matrix class.
 *
 * Used to represent a square 2D matrix of fixed size (size can be set during
 * construction).
 *
 * @tparam DATA_TYPE Data type of the matrix elements.
 */
template <typename DATA_TYPE> class Matrix {
private:
  /*! @brief Size of the square 2D matrix in 1 dimension. */
  const uint_fast32_t _size;

  /*! @brief Flat array containing the matrix elements. */
  std::vector<DATA_TYPE> _array;

public:
  /**
   * @brief Constructor.
   *
   * @param size Size of the square 2D matrix in 1 dimension.
   */
  inline Matrix(const uint_fast32_t size) : _size(size) {
    _array.resize(size * size);
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
    return _array[i * _size + j];
  }

  /**
   * @brief Invert the matrix using an auxiliary PLU decomposition.
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
   *    U^{-1}_{ij} = -\frac{1}{U_{jj}} \sum_{k=i}{j} U^{-1}_{ik} R_{kj}.
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
   */
  inline void plu_inverse() {

    // allocate the pivot array
    std::vector<uint_fast32_t> pivot_array(_size, 0);
    // alias the object using the label A (to be consistent with the notation
    // in the function documentation)
    Matrix &A = *this;

    // step 1: PLU decomposition of the original matrix
    for (uint_fast32_t i = 0; i < _size; ++i) {
      // find the next pivot
      uint_fast32_t imax = i;
      // note that we always use a double to store the element absolute value,
      // independent of the element type (std::abs() returns a double when
      // applied to a std::complex<double>)
      double Smax = std::abs(A(imax, i));
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        const double this_Smax = std::abs(A(j, i));
        if (this_Smax > Smax) {
          Smax = this_Smax;
          imax = j;
        }
      }
      // store the pivot in the pivot array
      pivot_array[i] = imax;
      // swap rows if necessary
      if (i != imax) {
        for (uint_fast32_t j = 0; j < _size; ++j) {
          // we need a temporary variable to now overwrite the original A_ij
          const DATA_TYPE temp = A(i, j);
          A(i, j) = A(imax, j);
          A(imax, j) = temp;
        }
      }
      // compute the inverse of A_ii to save on divisions
      const DATA_TYPE Aii_inv = 1. / A(i, i);
      // now multiply all elements below row i in column i with this value
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        A(j, i) *= Aii_inv;
      }
      // eliminate lower triangular elements in column i from the future U
      // matrix
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        // since we store the L matrix in the lower diagonal part of the matrix,
        // we only update columns k > i
        for (uint_fast32_t k = i + 1; k < _size; ++k) {
          A(j, k) -= A(j, i) * A(i, k);
        }
      }
    }

    // step 2: inversion of the U part of A
    for (uint_fast32_t i = 0; i < _size; ++i) {
      // compute the diagonal element of the inverse matrix for row i
      A(i, i) = 1. / A(i, i);
      // now use forward substitution to compute the other elements in this row
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        // we need a temporary variable, as A_ij also features in the loop below
        DATA_TYPE Aij = 0.;
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
    std::vector<DATA_TYPE> work(_size);
    // note that we need to iterate backwards for this algorithm
    for (uint_fast32_t jp1 = _size; jp1 > 0; --jp1) {
      // unsigned counters overflow when you decrement them below 0, so we
      // need a second variable to handle j = 0 in the loop condition
      const uint_fast32_t j = jp1 - 1;
      // initialise the new column with the values in the old column above the
      // diagonal
      for (uint_fast32_t i = 0; i < jp1; ++i) {
        work[i] = A(i, j);
      }
      // zero out the other elements in the new column
      for (uint_fast32_t i = jp1; i < _size; ++i) {
        work[i] = 0.;
      }
      // now apply the summation to eliminate unknowns in the column
      for (uint_fast32_t k = jp1; k < _size; ++k) {
        for (uint_fast32_t i = 0; i < _size; ++i) {
          work[i] -= A(i, k) * A(k, j);
        }
      }
      // we are done with the old column, overwrite it with the new column
      for (uint_fast32_t i = 0; i < _size; ++i) {
        A(i, j) = work[i];
      }
    }

    // step 4: exchange columns based on the pivot array
    // again, this is a backward iteration
    // note that we swap columns, since we multiply the inverse permutation
    // matrix on the other side
    for (uint_fast32_t jp1 = _size; jp1 > 0; --jp1) {
      // again, we need a second variable to handle unsigned integer overflow
      const uint_fast32_t j = jp1 - 1;
      const uint_fast32_t jp = pivot_array[j];
      // only swap columns if the pivot index is different
      if (jp != j) {
        for (uint_fast32_t i = 0; i < _size; ++i) {
          const DATA_TYPE temp = A(i, j);
          A(i, j) = A(i, jp);
          A(i, jp) = temp;
        }
      }
    }
  }
};
