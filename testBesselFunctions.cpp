#include "BesselFunctions.hpp"
#include <fstream>

int main(int argc, char **argv) {

  std::ofstream ofile("test_bessel.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const double x = 0.05 + 0.1 * i;
    double yr[10], yi[10], dyr[10], dyi[10];
    BesselFunctions::spherical_j_dj_complex_array(10, x, x, yr, yi, dyr, dyi);
    ofile << x << "\t" << yr[0] << "\t" << yi[0] << "\t" << dyr[0] << "\t"
          << dyi[0] << "\n";
  }

  return 0;
}
