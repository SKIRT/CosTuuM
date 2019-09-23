/**
 * @file CTM.cpp
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
#include "SpecialFunctions.hpp"
#include "TMatrix.hpp"

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

  // make sure 'rat' contains the right ratio if it is not 1
  if (abs(ratio_of_radii - 1.) > 1.e-8) {
    ratio_of_radii =
        SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
            axis_ratio);
  }

  // R_V is the equivalent sphere radius
  const float_type R_V = ratio_of_radii * axi;
  // we need a reasonable initial guess for the order of the spherical harmonics
  // the below expression provides this guess
  const float_type xev = 2. * M_PI * R_V / wavelength;
  uint_fast32_t nmax = static_cast<uint_fast32_t>(
      std::max(float_type(4.), xev + 4.05 * cbrt(xev)));

  // we need to find a maximum expansion order and number of Gauss-Legendre
  // quadrature points that leads to a converged T-matrix
  // To this end, we will start with a reasonable guess for the order, and then
  // keep track of how the accuracy of the scattering and extinction factors
  // changes as we first increase the order and then the number of quadrature
  // points
  // To simplify the calculation, we only compute the m=0 degree terms in the
  // T-matrix during the convergence loops
  TMatrix *active_Tmatrix = nullptr;

  // loop control variables
  float_type old_qext = 0.;
  float_type old_qsca = 0.;
  float_type dext = 1.;
  float_type dsca = 1.;
  while (nmax < maximum_order && (dext > tolerance || dsca > tolerance)) {
    // initially we assume that the number of quadrature points is a fixed
    // multiple of the order
    const uint_fast32_t ngauss = ndgs * nmax;

    // delete the old matrix (if it exists) and create a new one
    delete active_Tmatrix;
    active_Tmatrix = new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
    TMatrix &T = *active_Tmatrix;

    // calculate the scattering and extinction factors for this iteration
    float_type qsca = 0.;
    float_type qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const float_type dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    // compute the relative difference w.r.t. the previous values
    dsca = abs((old_qsca - qsca) / qsca);
    dext = abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    // some (temporary) diagnostic output
    ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

    // increase the order
    ++nmax;
  }

  // check if we managed to converge
  if (nmax == maximum_order) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in last iteration
    --nmax;
    ctm_warning("Converged for nmax = %" PRIuFAST32, nmax);
  }

  // now start increasing the number of quadrature points
  uint_fast32_t ngauss = nmax * ndgs + 1;
  dext = tolerance + 1.;
  dsca = tolerance + 1.;
  while (ngauss < maximum_ngauss && (dext > tolerance || dsca > tolerance)) {

    // delete the old matrix and create a new one
    delete active_Tmatrix;
    active_Tmatrix = new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
    TMatrix &T = *active_Tmatrix;

    // calculate the scattering and extinction factors
    float_type qsca = 0.;
    float_type qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const float_type dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    // compute the relative difference w.r.t. the old values
    dsca = abs((old_qsca - qsca) / qsca);
    dext = abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    // some diagnostic output
    ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

    // increase the number of quadrature points
    ++ngauss;
  }

  // check if we reached convergence
  if (ngauss == maximum_ngauss) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in final iteration
    --ngauss;
    ctm_warning("Converged for ngauss = %" PRIuFAST32, ngauss);
  }

  // ok: we found the right order and number of quadrature points
  // we now need to compute the additional elements for m=/=0
  active_Tmatrix->compute_additional_elements();

  // compute the actual scattering and extinction factors using the full matrix
  TMatrix &T = *active_Tmatrix;
  old_qsca = 0.;
  old_qext = 0.;
  for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
    for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
      for (int_fast32_t m1 = -n1; m1 < static_cast<int_fast32_t>(n1 + 1);
           ++m1) {
        for (int_fast32_t m2 = -n2; m2 < static_cast<int_fast32_t>(n2 + 1);
             ++m2) {
          old_qsca += std::norm(T(0, n1, m1, 0, n2, m2));
          old_qsca += std::norm(T(0, n1, m1, 1, n2, m2));
          old_qsca += std::norm(T(1, n1, m1, 0, n2, m2));
          old_qsca += std::norm(T(1, n1, m1, 1, n2, m2));
        }
      }
    }
  }
  for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
    for (int_fast32_t m1 = -n1; m1 < static_cast<int_fast32_t>(n1 + 1); ++m1) {
      old_qext += T(0, n1, m1, 0, n1, m1).real();
      old_qext += T(0, n1, m1, 1, n1, m1).real();
      old_qext += T(1, n1, m1, 0, n1, m1).real();
      old_qext += T(1, n1, m1, 1, n1, m1).real();
    }
  }
  // output the factors
  ctm_warning("qsca: %g", double(old_qsca));
  ctm_warning("qext: %g", double(old_qext));
  ctm_warning("walb: %g", double(-old_qsca / old_qext));

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

  Z.binary_dump(parser.get_value<std::string>("output"));

  // clean up
  delete active_Tmatrix;

  // done!
  return 0;
}
