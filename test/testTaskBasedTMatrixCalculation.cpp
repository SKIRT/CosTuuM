/**
 * @file testTaskBasedTMatrixCalculation.cpp
 *
 * @brief Unit test for the task based T-matrix calculation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "Error.hpp"
#include "GaussBasedResources.hpp"
#include "InteractionResource.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "TMatrixResource.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include "WignerDResources.hpp"

#include <cinttypes>
#include <fstream>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the task based T-matrix calculation.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const bool do_serial_test = true;
  const bool do_quicksched_test = false;

  std::ifstream ifile("test_tmatrixcalculator.txt");
  std::string line;
  // skip the first comment line
  getline(ifile, line);
  // now read the first data line
  getline(ifile, line);
  std::istringstream linestream(line);
  float_type axi, rat, lam, mrr, mri, eps, ddelt, alpha, beta, thet0, thet,
      phi0, phi, refqsca, refqext, refwalb, refZ[4][4];
  uint_fast32_t ndgs;
  linestream >> axi >> rat >> lam >> mrr >> mri >> eps >> ddelt >> ndgs >>
      alpha >> beta >> thet0 >> thet >> phi0 >> phi >> refqsca >> refqext >>
      refwalb >> refZ[0][0] >> refZ[0][1] >> refZ[0][2] >> refZ[0][3] >>
      refZ[1][0] >> refZ[1][1] >> refZ[1][2] >> refZ[1][3] >> refZ[2][0] >>
      refZ[2][1] >> refZ[2][2] >> refZ[2][3] >> refZ[3][0] >> refZ[3][1] >>
      refZ[3][2] >> refZ[3][3];

  const float_type particle_radius =
      UnitConverter::to_SI<QUANTITY_LENGTH>(double(axi), "micron");
  const float_type wavelength =
      UnitConverter::to_SI<QUANTITY_LENGTH>(double(lam), "micron");
  const float_type axis_ratio = eps;
  float_type ratio_of_radii;
  if (abs(rat - 1.) > 1.e-8) {
    ratio_of_radii =
        SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
            axis_ratio);
  } else {
    ratio_of_radii = rat;
  }
  // R_V is the equivalent sphere radius
  const float_type R_V = ratio_of_radii * particle_radius;
  const float_type xev = 2. * M_PI * R_V / wavelength;
  const float_type tolerance = ddelt * 0.1;
  const std::complex<float_type> refractive_index(mrr, mri);
  const uint_fast32_t maximum_order = 100;
  const uint_fast32_t minimum_order = static_cast<uint_fast32_t>(
      std::max(float_type(4.), xev + 4.05 * cbrt(xev)));
  const uint_fast32_t auxsize = 100;

  ctm_warning("Minimum order: %" PRIuFAST32, minimum_order);

  /// Serial version
  if (do_serial_test) {
    NBasedResources nfactors(maximum_order);
    nfactors.execute();

    std::vector<TMatrixAuxiliarySpace *> auxspace(auxsize, nullptr);
    for (uint_fast32_t i = 0; i < auxsize; ++i) {
      auxspace[i] = new TMatrixAuxiliarySpace(maximum_order);
    }

    std::vector<GaussBasedResources *> quadrature_points(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      quadrature_points[i] =
          new GaussBasedResources(ndgs * (minimum_order + i));
      quadrature_points[i]->execute();
    }

    std::vector<WignerDResources *> wignerdm0(maximum_order - minimum_order,
                                              nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      wignerdm0[i] =
          new WignerDResources(0, minimum_order + i, ndgs * (minimum_order + i),
                               *quadrature_points[i]);
      wignerdm0[i]->execute();
    }

    std::vector<ParticleGeometryResource *> geometries(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      geometries[i] = new ParticleGeometryResource(
          R_V, axis_ratio, ndgs * (minimum_order + i), *quadrature_points[i]);
      geometries[i]->execute();
    }

    std::vector<InteractionResource *> interactions(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      interactions[i] = new InteractionResource(
          wavelength, refractive_index, minimum_order + i,
          ndgs * (minimum_order + i), *geometries[i]);
      interactions[i]->execute();
    }

    TMatrixResource Tmatrix(maximum_order);
    std::vector<TMatrixM0Task *> m0tasks(maximum_order - minimum_order,
                                         nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      m0tasks[i] = new TMatrixM0Task(
          tolerance, minimum_order + i, ndgs * (minimum_order + i), nfactors,
          *quadrature_points[i], *geometries[i], *interactions[i],
          *wignerdm0[i], *auxspace[0], Tmatrix, Tmatrix.get_m_resource(0));
      m0tasks[i]->execute();
    }

    ctm_warning("nmax: %" PRIuFAST32, Tmatrix.get_nmax());
    ctm_warning("ngauss: %" PRIuFAST32, Tmatrix.get_ngauss());
    ctm_warning("Qsca: %g (%g)", double(Tmatrix.get_scattering_coefficient()),
                double(refqsca));
    ctm_warning("Qext: %g (%g)", double(Tmatrix.get_extinction_coefficient()),
                double(refqext));

    for (uint_fast32_t i = 0; i < auxsize; ++i) {
      delete auxspace[i];
    }
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      delete quadrature_points[i];
      delete wignerdm0[i];
      delete geometries[i];
      delete interactions[i];
      delete m0tasks[i];
    }
  }

  /// QuickSched version
  if (do_quicksched_test) {
    QuickSched quicksched(4);

    NBasedResources nfactors(maximum_order);
    quicksched.register_resource(nfactors);
    quicksched.register_task(nfactors);
    nfactors.link_resources(quicksched);

    std::vector<TMatrixAuxiliarySpace *> auxspace(auxsize, nullptr);
    for (uint_fast32_t i = 0; i < auxsize; ++i) {
      auxspace[i] = new TMatrixAuxiliarySpace(maximum_order);
      quicksched.register_resource(*auxspace[i]);
    }

    TMatrixResource Tmatrix(maximum_order);
    for (uint_fast32_t m = 0; m < maximum_order + 1; ++m) {
      quicksched.register_resource(Tmatrix.get_m_resource(m));
    }

    std::vector<GaussBasedResources *> quadrature_points(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      quadrature_points[i] =
          new GaussBasedResources(ndgs * (minimum_order + i));
      quicksched.register_resource(*quadrature_points[i]);
      quicksched.register_task(*quadrature_points[i]);
    }

    std::vector<WignerDResources *> wignerdm0(maximum_order - minimum_order,
                                              nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      wignerdm0[i] =
          new WignerDResources(0, minimum_order + i, ndgs * (minimum_order + i),
                               *quadrature_points[i]);
      quicksched.register_resource(*wignerdm0[i]);
      quicksched.register_task(*wignerdm0[i]);
      wignerdm0[i]->link_resources(quicksched);
      quicksched.link_tasks(*quadrature_points[i], *wignerdm0[i]);
    }

    std::vector<ParticleGeometryResource *> geometries(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      geometries[i] = new ParticleGeometryResource(
          R_V, axis_ratio, ndgs * (minimum_order + i), *quadrature_points[i]);
      quicksched.register_resource(*geometries[i]);
      quicksched.register_task(*geometries[i]);
      geometries[i]->link_resources(quicksched);
      quicksched.link_tasks(*quadrature_points[i], *geometries[i]);
    }

    std::vector<InteractionResource *> interactions(
        maximum_order - minimum_order, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      interactions[i] = new InteractionResource(
          wavelength, refractive_index, minimum_order + i,
          ndgs * (minimum_order + i), *geometries[i]);
      quicksched.register_resource(*interactions[i]);
      quicksched.register_task(*interactions[i]);
      interactions[i]->link_resources(quicksched);
      quicksched.link_tasks(*geometries[i], *interactions[i]);
    }

    std::vector<TMatrixM0Task *> m0tasks(maximum_order - minimum_order,
                                         nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      m0tasks[i] = new TMatrixM0Task(
          tolerance, minimum_order + i, ndgs * (minimum_order + i), nfactors,
          *quadrature_points[i], *geometries[i], *interactions[i],
          *wignerdm0[i], *auxspace[0], Tmatrix, Tmatrix.get_m_resource(0));
      quicksched.register_task(*m0tasks[i]);
      m0tasks[i]->link_resources(quicksched);
      quicksched.link_tasks(nfactors, *m0tasks[i]);
      quicksched.link_tasks(*wignerdm0[i], *m0tasks[i]);
      quicksched.link_tasks(*interactions[i], *m0tasks[i]);
      if (i > 0) {
        quicksched.link_tasks(*m0tasks[i - 1], *m0tasks[i]);
      }
    }

    quicksched.execute_tasks(4);

    ctm_warning("Qsca: %g", double(Tmatrix.get_scattering_coefficient()));
    ctm_warning("Qext: %g", double(Tmatrix.get_extinction_coefficient()));

    std::ofstream taskfile("test_tmatrix_tasks.txt");
    taskfile << "# thread\tstart\tend\ttype\n";
    quicksched.print_task(nfactors, taskfile);
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      quicksched.print_task(*quadrature_points[i], taskfile);
      quicksched.print_task(*wignerdm0[i], taskfile);
      quicksched.print_task(*geometries[i], taskfile);
      quicksched.print_task(*interactions[i], taskfile);
      quicksched.print_task(*m0tasks[i], taskfile);
    }
    std::ofstream typefile("test_tmatrix_types.txt");
    typefile << "# type\tlabel\n";
    quicksched.print_type_dict(typefile);

    // as so nicely stated throughout SWIFT: be clean
    for (uint_fast32_t i = 0; i < auxsize; ++i) {
      delete auxspace[i];
    }
    for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
      delete quadrature_points[i];
      delete wignerdm0[i];
      delete geometries[i];
      delete interactions[i];
      delete m0tasks[i];
    }
  }

  return 0;
}
