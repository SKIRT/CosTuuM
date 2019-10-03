/**
 * @file CosTuuM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "CommandLineParser.hpp"
#include "Configuration.hpp"
#include "ConfigurationInfo.hpp"
#include "Matrix.hpp"
#include "ParameterFile.hpp"
#include "TMatrixCalculator.hpp"

#include <cinttypes>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Main program entry point.
 *
 * For now, does the entire T-matrix calculation for the Mishchenko default
 * parameters and outputs the scattering and extinction factors.
 *
 * @param argc Number of command line arguments (currently ignored).
 * @param argv Command line arguments (currently ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  std::cout << "Program configuration:" << std::endl;
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::cout << it.get_key() << ": " << it.get_value() << "\n";
  }
  std::cout << std::endl;

  CommandLineParser parser("C++ TMatrix code");
  parser.add_option("params", 'p', "Name of the parameter file.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "", true);
  parser.add_option("output", 'o', "Name of the output dump file.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "", true);
  parser.parse_arguments(argc, argv);

  std::cout << "Command line options:" << std::endl;
  parser.print_contents(std::cout);
  std::cout << std::endl;

  ParameterFile params(parser.get_value<std::string>("params"));

  /// input parameters

  /// dust particle
  // size of the particle (in same units as the wavelength)
  const float_type axi = params.get_physical_value<QUANTITY_LENGTH>(
                             "DustParticle:size", "10. micron") *
                         1.e6;
  // ratio between the equal surface area sphere radius and equal volume sphere
  // radius (is recomputed if not equal to 1)
  float_type ratio_of_radii = 1.;
  if (!params.get_value<bool>("DustParticle:is equal volume sphere radius",
                              true)) {
    ratio_of_radii = 0.1;
  }
  // refractive index
  const std::complex<float_type> mr = params.get_value<std::complex<double>>(
      "DustParticle:refractive index", std::complex<double>(1.5, 0.02));
  // ratio of horizontal and vertical axis of the spheroidal particle
  const float_type axis_ratio =
      params.get_value<double>("DustParticle:axis ratio", 0.5);

  /// radiation
  // wavelength of incoming radiation (in same units as the particle size)
  const float_type wavelength =
      params.get_physical_value<QUANTITY_LENGTH>(
          "Radiation:incoming wavelength", "6.283185307 micron") *
      1.e6;

  /// calculation
  // tolerance for the calculation
  const float_type tolerance =
      params.get_value<double>("Calculation:tolerance", 1.e-4);
  // number of Gauss-Legendre points to use as a multiplier of the maximum
  // order of spherical harmonics used to decompose the electric field, during
  // the first loop of the algorithm
  const uint_fast32_t ndgs =
      params.get_value<uint_fast32_t>("Calculation:Gauss Legendre factor", 2);
  // maximum number of iterations during the first loop
  const uint_fast32_t maximum_order =
      params.get_value<uint_fast32_t>("Calculation:maximum order", 200);
  // maximum number of iterations during the second loop
  const uint_fast32_t maximum_ngauss = params.get_value<uint_fast32_t>(
      "Calculation:maximum number of Gauss Legendre points", 500);

  /// scattering event
  // first and second Euler angles describing the orientation of the particle
  // in the fixed laboratory reference frame
  const float_type alpha = params.get_physical_value<QUANTITY_ANGLE>(
      "DustParticle:alpha angle", "145. degrees");
  const float_type beta = params.get_physical_value<QUANTITY_ANGLE>(
      "DustParticle:beta angle", "52. degrees");

  // zenith and azimuth angle of the incoming photon
  const float_type theta_in = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:incoming theta", "56. degrees");
  const float_type phi_in = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:incoming phi", "114. degrees");

  // zenith and azimuth angle of the scattered photon
  const float_type theta_out = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:outgoing theta", "65. degrees");
  const float_type phi_out = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:outgoing phi", "128. degrees");

  // write used parameters to file
  {
    std::ofstream pfile("parameters-usedvalues.param");
    params.print_contents(pfile);
    pfile.close();

    std::cout << "Input parameters:" << std::endl;
    params.print_contents(std::cout);
    std::cout << std::endl;
  }

  TMatrix *active_Tmatrix = TMatrixCalculator::calculate_TMatrix(
      ratio_of_radii, axis_ratio, axi, wavelength, maximum_order, tolerance,
      ndgs, mr, maximum_ngauss);

  /// compute a scattering event using the T-matrix

  Matrix<float_type> Z = active_Tmatrix->get_scattering_matrix(
      alpha, beta, theta_in, phi_in, theta_out, phi_out);

  ctm_warning("Z[0,:] = %g %g %g %g", double(Z(0, 0)), double(Z(0, 1)),
              double(Z(0, 2)), double(Z(0, 3)));
  ctm_warning("Z[1,:] = %g %g %g %g", double(Z(1, 0)), double(Z(1, 1)),
              double(Z(1, 2)), double(Z(1, 3)));
  ctm_warning("Z[2,:] = %g %g %g %g", double(Z(2, 0)), double(Z(2, 1)),
              double(Z(2, 2)), double(Z(2, 3)));
  ctm_warning("Z[3,:] = %g %g %g %g", double(Z(3, 0)), double(Z(3, 1)),
              double(Z(3, 2)), double(Z(3, 3)));

  Matrix<float_type> K =
      active_Tmatrix->get_extinction_matrix(alpha, beta, theta_in, phi_in);

  ctm_warning("K[0,:] = %g %g %g %g", double(K(0, 0)), double(K(0, 1)),
              double(K(0, 2)), double(K(0, 3)));
  ctm_warning("K[1,:] = %g %g %g %g", double(K(1, 0)), double(K(1, 1)),
              double(K(1, 2)), double(K(1, 3)));
  ctm_warning("K[2,:] = %g %g %g %g", double(K(2, 0)), double(K(2, 1)),
              double(K(2, 2)), double(K(2, 3)));
  ctm_warning("K[3,:] = %g %g %g %g", double(K(3, 0)), double(K(3, 1)),
              double(K(3, 2)), double(K(3, 3)));

  Z.binary_dump(parser.get_value<std::string>("output"));

  // clean up
  delete active_Tmatrix;

  // done!
  return 0;
}
