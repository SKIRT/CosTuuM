#include "Matrix.hpp"
#include <iostream>

int main(int argc, char **argv) {

  Matrix A(3);
  A(0, 0) = 1.;
  A(0, 0).imag(1.);
  A(0, 1) = 2.;
  A(0, 1).imag(1.);
  A(0, 2) = 3.;
  A(0, 2).imag(1.);
  A(1, 0) = 4.;
  A(1, 0).imag(2.);
  A(1, 1) = 5.;
  A(1, 1).imag(1.);
  A(1, 2) = 6.;
  A(1, 2).imag(1.);
  A(2, 0) = 8.;
  A(2, 0).imag(3.);
  A(2, 1) = 9.;
  A(2, 1).imag(1.);
  A(2, 2).imag(-2.);

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  A.plu_inverse();

  for (uint_fast32_t i = 0; i < 3; ++i) {
    for (uint_fast32_t j = 0; j < 3; ++j) {
      std::cout << A(i, j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
