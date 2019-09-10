#include <cinttypes>
#include <complex>
#include <vector>

class Matrix {
private:
  const uint_fast32_t _size;

  std::vector<std::complex<double>> _array;

public:
  inline Matrix(const uint_fast32_t size) : _size(size) {
    _array.resize(size * size);
  }

  inline std::complex<double> &operator()(const uint_fast32_t i,
                                          const uint_fast32_t j) {
    return _array[i * _size + j];
  }

  inline void plu_inverse() {

    std::vector<uint_fast32_t> pivot_array(_size, 0);
    Matrix &A = *this;

    // step 1: PLU decomposition of the original matrix
    for (uint_fast32_t i = 0; i < _size; ++i) {
      uint_fast32_t imax = i;
      double Smax = std::abs(A(imax, i));
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        const double this_Smax = std::abs(A(j, i));
        if (this_Smax > Smax) {
          Smax = this_Smax;
          imax = j;
        }
      }
      pivot_array[i] = imax;
      if (i != imax) {
        for (uint_fast32_t j = 0; j < _size; ++j) {
          const std::complex<double> temp = A(i, j);
          A(i, j) = A(imax, j);
          A(imax, j) = temp;
        }
      }
      const std::complex<double> Aii_inv = 1. / A(i, i);
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        A(j, i) *= Aii_inv;
      }
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        for (uint_fast32_t k = i + 1; k < _size; ++k) {
          A(j, k) -= A(j, i) * A(i, k);
        }
      }
    }

    // step 2: inversion of the U part of A
    for (uint_fast32_t i = 0; i < _size; ++i) {
      A(i, i) = 1. / A(i, i);
      for (uint_fast32_t j = i + 1; j < _size; ++j) {
        std::complex<double> Aij = 0.;
        for (uint_fast32_t k = i; k < j; ++k) {
          Aij -= A(i, k) * A(k, j);
        }
        A(i, j) = Aij / A(j, j);
      }
    }

    // step 3: solve inv(A)*L = inv(U)
    std::vector<std::complex<double>> work(_size);
    for (uint_fast32_t jp1 = _size; jp1 > 0; --jp1) {
      const uint_fast32_t j = jp1 - 1;
      for (uint_fast32_t i = 0; i < jp1; ++i) {
        work[i] = A(i, j);
      }
      for (uint_fast32_t i = jp1; i < _size; ++i) {
        work[i] = 0.;
      }
      for (uint_fast32_t k = jp1; k < _size; ++k) {
        for (uint_fast32_t i = 0; i < _size; ++i) {
          work[i] -= A(i, k) * A(k, j);
        }
      }
      for (uint_fast32_t i = 0; i < _size; ++i) {
        A(i, j) = work[i];
      }
    }

    // step 4: exchange columns based on the pivot matrix
    for (uint_fast32_t jp1 = _size; jp1 > 0; --jp1) {
      const uint_fast32_t j = jp1 - 1;
      const uint_fast32_t jp = pivot_array[j];
      if (jp != j) {
        for (uint_fast32_t i = 0; i < _size; ++i) {
          const std::complex<double> temp = A(i, j);
          A(i, j) = A(i, jp);
          A(i, jp) = temp;
        }
      }
    }
  }
};
