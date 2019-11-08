/**
 * @file testTaskBasedTMatrixCalculation.cpp
 *
 * @brief Unit test for the task based T-matrix calculation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "Configuration.hpp"
#include "ConvergedSizeResources.hpp"
#include "Error.hpp"
#include "ExtinctionMatrixResource.hpp"
#include "GaussBasedResources.hpp"
#include "InteractionResource.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "ScatteringMatrixResource.hpp"
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
 * @brief Delete all pointers in the given pointer vector.
 *
 * @param vec Pointer vector.
 */
template <typename T> inline void clear_vector(std::vector<T *> &vec) {
  for (uint_fast32_t i = 0; i < vec.size(); ++i) {
    delete vec[i];
  }
}

/**
 * @brief Print all tasks in the given vector to the given file.
 *
 * @param vec Vector containing task pointers.
 * @param quicksched QuickSched library wrapper.
 * @param ofile Output file.
 */
template <typename T>
inline void print_vector(std::vector<T *> &vec, QuickSched &quicksched,
                         std::ofstream &ofile) {
  for (uint_fast32_t i = 0; i < vec.size(); ++i) {
    quicksched.print_task(*vec[i], ofile);
  }
}

/**
 * @brief Unit test for the task based T-matrix calculation.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const bool do_serial_test = true;
  const bool do_quicksched_test = true;
  const bool do_full_benchmark_test = true;

  if (do_serial_test || do_quicksched_test) {
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

    axi = UnitConverter::to_SI<QUANTITY_LENGTH>(double(axi), "micron");
    lam = UnitConverter::to_SI<QUANTITY_LENGTH>(double(lam), "micron");
    alpha = UnitConverter::to_SI<QUANTITY_ANGLE>(double(alpha), "degrees");
    beta = UnitConverter::to_SI<QUANTITY_ANGLE>(double(beta), "degrees");
    thet0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet0), "degrees");
    thet = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet), "degrees");
    phi0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi0), "degrees");
    phi = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi), "degrees");
    // convert the reference results to SI units
    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 4; ++j) {
        refZ[i][j] = UnitConverter::to_SI<QUANTITY_SURFACE_AREA>(
            double(refZ[i][j]), "micron^2");
      }
    }

    const float_type particle_radius = axi;
    const float_type wavelength = lam;
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

      ConvergedSizeResources converged_size;

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

      std::vector<WignerDm0Resources *> wignerdm0(maximum_order - minimum_order,
                                                  nullptr);
      for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
        wignerdm0[i] = new WignerDm0Resources(minimum_order + i,
                                              ndgs * (minimum_order + i),
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
            *wignerdm0[i], *auxspace[0], Tmatrix, converged_size,
            Tmatrix.get_m_resource(0));
        m0tasks[i]->execute();
      }

      ctm_warning("nmax: %" PRIuFAST32, Tmatrix.get_nmax());
      ctm_warning("ngauss: %" PRIuFAST32, Tmatrix.get_ngauss());
      ctm_warning("Qsca: %g (%g)", double(Tmatrix.get_scattering_coefficient()),
                  double(refqsca));
      ctm_warning("Qext: %g (%g)", double(Tmatrix.get_extinction_coefficient()),
                  double(refqext));

      std::vector<WignerDmn0Resources *> more_wigner(maximum_order, nullptr);
      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        more_wigner[i] = new WignerDmn0Resources(
            1 + i, maximum_order, ndgs * maximum_order, converged_size);
        more_wigner[i]->execute();
      }

      std::vector<TMatrixMAllTask *> malltask(maximum_order, nullptr);
      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        malltask[i] = new TMatrixMAllTask(
            1 + i, nfactors, *more_wigner[i], converged_size,
            *auxspace[i % auxsize], Tmatrix, Tmatrix.get_m_resource(1 + i));
        malltask[i]->execute();
      }

      TMatrixQTask qtask(Tmatrix, Tmatrix.get_m_resource(0));
      qtask.execute();

      const float_type qext = Tmatrix.get_extinction_coefficient();
      const float_type qsca = Tmatrix.get_scattering_coefficient();
      const float_type walb = -qsca / qext;
      ctm_warning("Qsca: %g (%g)", double(qsca), double(refqsca));
      ctm_warning("Qext: %g (%g)", double(qext), double(refqext));

      assert_values_equal_rel(double(qext), double(refqext), 1.e-5);
      assert_values_equal_rel(double(qsca), double(refqsca), 1.e-5);
      assert_values_equal_rel(double(walb), double(refwalb), 1.e-5);

      ScatteringMatrixResource Zmatrix(alpha, beta, thet0, phi0, thet, phi,
                                       Tmatrix, maximum_order);
      Zmatrix.execute();

      ctm_warning("Z[0,:]: %g %g %g %g", double(Zmatrix(0, 0)),
                  double(Zmatrix(0, 1)), double(Zmatrix(0, 2)),
                  double(Zmatrix(0, 3)));
      ctm_warning("Zref[0,:]: %g %g %g %g", double(refZ[0][0]),
                  double(refZ[0][1]), double(refZ[0][2]), double(refZ[0][3]));
      ctm_warning("Z[1,:]: %g %g %g %g", double(Zmatrix(1, 0)),
                  double(Zmatrix(1, 1)), double(Zmatrix(1, 2)),
                  double(Zmatrix(1, 3)));
      ctm_warning("Zref[1,:]: %g %g %g %g", double(refZ[1][0]),
                  double(refZ[1][1]), double(refZ[1][2]), double(refZ[1][3]));
      ctm_warning("Z[2,:]: %g %g %g %g", double(Zmatrix(2, 0)),
                  double(Zmatrix(2, 1)), double(Zmatrix(2, 2)),
                  double(Zmatrix(2, 3)));
      ctm_warning("Zref[2,:]: %g %g %g %g", double(refZ[2][0]),
                  double(refZ[2][1]), double(refZ[2][2]), double(refZ[2][3]));
      ctm_warning("Z[3,:]: %g %g %g %g", double(Zmatrix(3, 0)),
                  double(Zmatrix(3, 1)), double(Zmatrix(3, 2)),
                  double(Zmatrix(3, 3)));
      ctm_warning("Zref[3,:]: %g %g %g %g", double(refZ[3][0]),
                  double(refZ[3][1]), double(refZ[3][2]), double(refZ[3][3]));

      // compare the result with the reference
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          assert_values_equal_rel(double(Zmatrix(i, j)), double(refZ[i][j]),
                                  2.e-2);
        }
      }

      clear_vector(auxspace);
      clear_vector(quadrature_points);
      clear_vector(wignerdm0);
      clear_vector(geometries);
      clear_vector(interactions);
      clear_vector(m0tasks);
      clear_vector(more_wigner);
      clear_vector(malltask);
    }

    /// QuickSched version
    if (do_quicksched_test) {
      QuickSched quicksched(4, true, "quicksched.log");

      NBasedResources nfactors(maximum_order);
      quicksched.register_resource(nfactors);
      quicksched.register_task(nfactors);
      nfactors.link_resources(quicksched);

      ConvergedSizeResources converged_size;
      quicksched.register_resource(converged_size);

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

      std::vector<WignerDm0Resources *> wignerdm0(maximum_order - minimum_order,
                                                  nullptr);
      for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
        wignerdm0[i] = new WignerDm0Resources(minimum_order + i,
                                              ndgs * (minimum_order + i),
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
            *wignerdm0[i], *auxspace[0], Tmatrix, converged_size,
            Tmatrix.get_m_resource(0));
        quicksched.register_task(*m0tasks[i]);
        m0tasks[i]->link_resources(quicksched);
        quicksched.link_tasks(nfactors, *m0tasks[i]);
        quicksched.link_tasks(*wignerdm0[i], *m0tasks[i]);
        quicksched.link_tasks(*interactions[i], *m0tasks[i]);
        if (i > 0) {
          quicksched.link_tasks(*m0tasks[i - 1], *m0tasks[i]);
        }
      }

      std::vector<WignerDmn0Resources *> more_wigner(maximum_order, nullptr);
      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        more_wigner[i] = new WignerDmn0Resources(
            1 + i, maximum_order, ndgs * maximum_order, converged_size);
        quicksched.register_resource(*more_wigner[i]);
        quicksched.register_task(*more_wigner[i]);
        more_wigner[i]->link_resources(quicksched);
        quicksched.link_tasks(*m0tasks.back(), *more_wigner[i]);
      }

      TMatrixQTask qtask(Tmatrix, Tmatrix.get_m_resource(0));
      quicksched.register_task(qtask);
      qtask.link_resources(quicksched);
      ExtinctionMatrixResource Kmatrix(0.3 * M_PI, 0., Tmatrix);
      quicksched.register_resource(Kmatrix);
      quicksched.register_task(Kmatrix);
      Kmatrix.link_resources(quicksched);
      ScatteringMatrixResource Zmatrix(alpha, beta, thet0, phi0, thet, phi,
                                       Tmatrix, maximum_order);
      quicksched.register_resource(Zmatrix);
      quicksched.register_task(Zmatrix);
      Zmatrix.link_resources(quicksched);

      std::vector<TMatrixMAllTask *> malltask(maximum_order, nullptr);
      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        malltask[i] = new TMatrixMAllTask(
            1 + i, nfactors, *more_wigner[i], converged_size,
            *auxspace[i % auxsize], Tmatrix, Tmatrix.get_m_resource(1 + i));
        quicksched.register_task(*malltask[i]);
        malltask[i]->link_resources(quicksched);
        quicksched.link_tasks(*more_wigner[i], *malltask[i]);

        quicksched.link_tasks(*malltask[i], qtask);
        quicksched.link_tasks(*malltask[i], Kmatrix);
        quicksched.link_tasks(*malltask[i], Zmatrix);
      }

      quicksched.execute_tasks(4);

      ctm_warning("nmax: %" PRIuFAST32, Tmatrix.get_nmax());
      ctm_warning("ngauss: %" PRIuFAST32, Tmatrix.get_ngauss());
      ctm_warning("Qsca: %g", double(Tmatrix.get_scattering_coefficient()));
      ctm_warning("Qext: %g", double(Tmatrix.get_extinction_coefficient()));

      const float_type qext = Tmatrix.get_extinction_coefficient();
      const float_type qsca = Tmatrix.get_scattering_coefficient();
      const float_type walb = -qsca / qext;

      assert_values_equal_rel(double(qext), double(refqext), 1.e-5);
      assert_values_equal_rel(double(qsca), double(refqsca), 1.e-5);
      assert_values_equal_rel(double(walb), double(refwalb), 1.e-5);

      ctm_warning("Qsca: %g (%g)", double(qsca), double(refqsca));
      ctm_warning("Qext: %g (%g)", double(qext), double(refqext));

      ctm_warning("K[0,:]: %g %g %g %g", double(Kmatrix(0, 0)),
                  double(Kmatrix(0, 1)), double(Kmatrix(0, 2)),
                  double(Kmatrix(0, 3)));
      ctm_warning("K[1,:]: %g %g %g %g", double(Kmatrix(1, 0)),
                  double(Kmatrix(1, 1)), double(Kmatrix(1, 2)),
                  double(Kmatrix(1, 3)));
      ctm_warning("K[2,:]: %g %g %g %g", double(Kmatrix(2, 0)),
                  double(Kmatrix(2, 1)), double(Kmatrix(2, 2)),
                  double(Kmatrix(2, 3)));
      ctm_warning("K[3,:]: %g %g %g %g", double(Kmatrix(3, 0)),
                  double(Kmatrix(3, 1)), double(Kmatrix(3, 2)),
                  double(Kmatrix(3, 3)));

      ctm_warning("Z[0,:]: %g %g %g %g", double(Zmatrix(0, 0)),
                  double(Zmatrix(0, 1)), double(Zmatrix(0, 2)),
                  double(Zmatrix(0, 3)));
      ctm_warning("Zref[0,:]: %g %g %g %g", double(refZ[0][0]),
                  double(refZ[0][1]), double(refZ[0][2]), double(refZ[0][3]));
      ctm_warning("Z[1,:]: %g %g %g %g", double(Zmatrix(1, 0)),
                  double(Zmatrix(1, 1)), double(Zmatrix(1, 2)),
                  double(Zmatrix(1, 3)));
      ctm_warning("Zref[1,:]: %g %g %g %g", double(refZ[1][0]),
                  double(refZ[1][1]), double(refZ[1][2]), double(refZ[1][3]));
      ctm_warning("Z[2,:]: %g %g %g %g", double(Zmatrix(2, 0)),
                  double(Zmatrix(2, 1)), double(Zmatrix(2, 2)),
                  double(Zmatrix(2, 3)));
      ctm_warning("Zref[2,:]: %g %g %g %g", double(refZ[2][0]),
                  double(refZ[2][1]), double(refZ[2][2]), double(refZ[2][3]));
      ctm_warning("Z[3,:]: %g %g %g %g", double(Zmatrix(3, 0)),
                  double(Zmatrix(3, 1)), double(Zmatrix(3, 2)),
                  double(Zmatrix(3, 3)));
      ctm_warning("Zref[3,:]: %g %g %g %g", double(refZ[3][0]),
                  double(refZ[3][1]), double(refZ[3][2]), double(refZ[3][3]));

      // compare the result with the reference
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          assert_values_equal_rel(double(Zmatrix(i, j)), double(refZ[i][j]),
                                  2.e-2);
        }
      }

      std::ofstream taskfile("test_tmatrix_tasks.txt");
      taskfile << "# thread\tstart\tend\ttype\ttask id\n";
      quicksched.print_task(nfactors, taskfile);
      for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
        quicksched.print_task(*quadrature_points[i], taskfile);
        quicksched.print_task(*wignerdm0[i], taskfile);
        quicksched.print_task(*geometries[i], taskfile);
        quicksched.print_task(*interactions[i], taskfile);
        quicksched.print_task(*m0tasks[i], taskfile);
      }
      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        quicksched.print_task(*more_wigner[i], taskfile);
        quicksched.print_task(*malltask[i], taskfile);
      }
      quicksched.print_task(qtask, taskfile);
      quicksched.print_task(Kmatrix, taskfile);
      quicksched.print_task(Zmatrix, taskfile);
      std::ofstream typefile("test_tmatrix_types.txt");
      typefile << "# type\tlabel\n";
      quicksched.print_type_dict(typefile);

      // as so nicely stated throughout SWIFT: be clean
      clear_vector(auxspace);
      clear_vector(quadrature_points);
      clear_vector(wignerdm0);
      clear_vector(geometries);
      clear_vector(interactions);
      clear_vector(m0tasks);
      clear_vector(more_wigner);
      clear_vector(malltask);
    }
  }

  // full comparison between the task based T-matrix and the Mishchenko T-matrix
  // code
  if (do_full_benchmark_test) {

    const uint_fast32_t auxsize = 100;
    const uint_fast32_t maximum_order = 40;
    const uint_fast32_t ndgs = 2;

    QuickSched quicksched(4, true, "test_Tmatrix_benchmark.log");

    NBasedResources nfactors(maximum_order);
    quicksched.register_resource(nfactors);
    quicksched.register_task(nfactors);
    nfactors.link_resources(quicksched);

    std::vector<TMatrixAuxiliarySpace *> auxspace(auxsize, nullptr);
    for (uint_fast32_t i = 0; i < auxsize; ++i) {
      auxspace[i] = new TMatrixAuxiliarySpace(maximum_order);
      quicksched.register_resource(*auxspace[i]);
    }

    std::vector<GaussBasedResources *> quadrature_points(maximum_order - 1,
                                                         nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - 1; ++i) {
      quadrature_points[i] = new GaussBasedResources(ndgs * (1 + i));
      quicksched.register_resource(*quadrature_points[i]);
      quicksched.register_task(*quadrature_points[i]);
    }

    std::vector<WignerDm0Resources *> wignerdm0(maximum_order - 1, nullptr);
    for (uint_fast32_t i = 0; i < maximum_order - 1; ++i) {
      wignerdm0[i] =
          new WignerDm0Resources(1 + i, ndgs * (1 + i), *quadrature_points[i]);
      quicksched.register_resource(*wignerdm0[i]);
      quicksched.register_task(*wignerdm0[i]);
      wignerdm0[i]->link_resources(quicksched);
      quicksched.link_tasks(*quadrature_points[i], *wignerdm0[i]);
    }

    std::vector<ConvergedSizeResources *> converged_sizes;
    std::vector<TMatrixResource *> Tmatrices;
    std::vector<TMatrixQTask *> Qtasks;
    std::vector<ScatteringMatrixResource *> Zmatrices;
    std::vector<ParticleGeometryResource *> geometries;
    std::vector<InteractionResource *> interactions;
    std::vector<TMatrixM0Task *> m0tasks;
    std::vector<WignerDmn0Resources *> wignerdmalls;
    std::vector<TMatrixMAllTask *> malltasks;
    std::vector<float_type> refQscavec;
    std::vector<float_type> refQextvec;
    std::vector<float_type> refwalbvec;
    std::vector<Matrix<float_type>> refZvec;
    uint_fast32_t auxcount = 0;
    std::ifstream ifile("test_tmatrixcalculator.txt");
    std::string line;
    // skip the first comment line
    getline(ifile, line);
    // now read the first data line
    while (getline(ifile, line)) {
      std::istringstream linestream(line);
      float_type axi, rat, lam, mrr, mri, eps, ddelt, alpha, beta, thet0, thet,
          phi0, phi, refqsca, refqext, refwalb;
      Matrix<float_type> refZ(4, 4);
      uint_fast32_t ndgs;
      linestream >> axi >> rat >> lam >> mrr >> mri >> eps >> ddelt >> ndgs >>
          alpha >> beta >> thet0 >> thet >> phi0 >> phi >> refqsca >> refqext >>
          refwalb >> refZ(0, 0) >> refZ(0, 1) >> refZ(0, 2) >> refZ(0, 3) >>
          refZ(1, 0) >> refZ(1, 1) >> refZ(1, 2) >> refZ(1, 3) >> refZ(2, 0) >>
          refZ(2, 1) >> refZ(2, 2) >> refZ(2, 3) >> refZ(3, 0) >> refZ(3, 1) >>
          refZ(3, 2) >> refZ(3, 3);

      axi = UnitConverter::to_SI<QUANTITY_LENGTH>(double(axi), "micron");
      lam = UnitConverter::to_SI<QUANTITY_LENGTH>(double(lam), "micron");
      alpha = UnitConverter::to_SI<QUANTITY_ANGLE>(double(alpha), "degrees");
      beta = UnitConverter::to_SI<QUANTITY_ANGLE>(double(beta), "degrees");
      thet0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet0), "degrees");
      thet = UnitConverter::to_SI<QUANTITY_ANGLE>(double(thet), "degrees");
      phi0 = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi0), "degrees");
      phi = UnitConverter::to_SI<QUANTITY_ANGLE>(double(phi), "degrees");
      // convert the reference results to SI units
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          refZ(i, j) = UnitConverter::to_SI<QUANTITY_SURFACE_AREA>(
              double(refZ(i, j)), "micron^2");
        }
      }
      refZvec.push_back(refZ);
      refQscavec.push_back(refqsca);
      refQextvec.push_back(refqext);
      refwalbvec.push_back(refwalb);

      const float_type particle_radius = axi;
      const float_type wavelength = lam;
      const float_type axis_ratio = eps;
      float_type ratio_of_radii;
      if (abs(rat - 1.) > 1.e-8) {
        ratio_of_radii = SpecialFunctions::
            get_equal_volume_to_equal_surface_area_sphere_ratio(axis_ratio);
      } else {
        ratio_of_radii = rat;
      }
      // R_V is the equivalent sphere radius
      const float_type R_V = ratio_of_radii * particle_radius;
      const float_type xev = 2. * M_PI * R_V / wavelength;
      const float_type tolerance = ddelt * 0.1;
      const std::complex<float_type> refractive_index(mrr, mri);
      const uint_fast32_t minimum_order = static_cast<uint_fast32_t>(
          std::max(float_type(4.), xev + 4.05 * cbrt(xev)));

      converged_sizes.push_back(new ConvergedSizeResources());
      quicksched.register_resource(*converged_sizes.back());
      Tmatrices.push_back(new TMatrixResource(maximum_order));
      for (uint_fast32_t m = 0; m < maximum_order + 1; ++m) {
        quicksched.register_resource(Tmatrices.back()->get_m_resource(m));
      }

      for (uint_fast32_t i = 0; i < maximum_order - minimum_order; ++i) {
        geometries.push_back(new ParticleGeometryResource(
            R_V, axis_ratio, ndgs * (minimum_order + i),
            *quadrature_points[minimum_order + i - 1]));
        quicksched.register_resource(*geometries.back());
        quicksched.register_task(*geometries.back());
        geometries.back()->link_resources(quicksched);
        quicksched.link_tasks(*quadrature_points[minimum_order + i - 1],
                              *geometries.back());

        interactions.push_back(new InteractionResource(
            wavelength, refractive_index, minimum_order + i,
            ndgs * (minimum_order + i), *geometries.back()));
        quicksched.register_resource(*interactions.back());
        quicksched.register_task(*interactions.back());
        interactions.back()->link_resources(quicksched);
        quicksched.link_tasks(*geometries.back(), *interactions.back());

        m0tasks.push_back(new TMatrixM0Task(
            tolerance, minimum_order + i, ndgs * (minimum_order + i), nfactors,
            *quadrature_points[minimum_order + i - 1], *geometries.back(),
            *interactions.back(), *wignerdm0[minimum_order + i - 1],
            *auxspace[auxcount], *Tmatrices.back(), *converged_sizes.back(),
            Tmatrices.back()->get_m_resource(0)));
        quicksched.register_task(*m0tasks.back());
        m0tasks.back()->link_resources(quicksched);
        quicksched.link_tasks(nfactors, *m0tasks.back());
        quicksched.link_tasks(*wignerdm0[minimum_order + i - 1],
                              *m0tasks.back());
        quicksched.link_tasks(*interactions.back(), *m0tasks.back());
        if (i > 0) {
          quicksched.link_tasks(*m0tasks[m0tasks.size() - 2], *m0tasks.back());
        }
      }

      Qtasks.push_back(new TMatrixQTask(*Tmatrices.back(),
                                        Tmatrices.back()->get_m_resource(0)));
      quicksched.register_task(*Qtasks.back());
      Qtasks.back()->link_resources(quicksched);
      Zmatrices.push_back(
          new ScatteringMatrixResource(alpha, beta, thet0, phi0, thet, phi,
                                       *Tmatrices.back(), maximum_order));
      quicksched.register_resource(*Zmatrices.back());
      quicksched.register_task(*Zmatrices.back());
      Zmatrices.back()->link_resources(quicksched);

      for (uint_fast32_t i = 0; i < maximum_order; ++i) {
        wignerdmalls.push_back(
            new WignerDmn0Resources(1 + i, maximum_order, ndgs * maximum_order,
                                    *converged_sizes.back()));
        quicksched.register_resource(*wignerdmalls.back());
        quicksched.register_task(*wignerdmalls.back());
        wignerdmalls.back()->link_resources(quicksched);
        quicksched.link_tasks(*m0tasks.back(), *wignerdmalls.back());

        malltasks.push_back(new TMatrixMAllTask(
            1 + i, nfactors, *wignerdmalls.back(), *converged_sizes.back(),
            *auxspace[auxcount], *Tmatrices.back(),
            Tmatrices.back()->get_m_resource(1 + i)));
        quicksched.register_task(*malltasks.back());
        malltasks.back()->link_resources(quicksched);
        quicksched.link_tasks(*wignerdmalls.back(), *malltasks.back());

        quicksched.link_tasks(*malltasks.back(), *Qtasks.back());
        quicksched.link_tasks(*malltasks.back(), *Zmatrices.back());
      }

      auxcount = (auxcount + 1) % auxsize;
    }

    quicksched.execute_tasks(4);

    std::ofstream taskfile("test_tmatrix_benchmark_tasks.txt");
    taskfile << "# thread\tstart\tend\ttype\ttask id\n";
    quicksched.print_task(nfactors, taskfile);
    print_vector(quadrature_points, quicksched, taskfile);
    print_vector(wignerdm0, quicksched, taskfile);
    print_vector(geometries, quicksched, taskfile);
    print_vector(interactions, quicksched, taskfile);
    print_vector(m0tasks, quicksched, taskfile);
    print_vector(wignerdmalls, quicksched, taskfile);
    print_vector(malltasks, quicksched, taskfile);
    print_vector(Qtasks, quicksched, taskfile);
    print_vector(Zmatrices, quicksched, taskfile);
    std::ofstream typefile("test_tmatrix_benchmark_types.txt");
    typefile << "# type\tlabel\n";
    quicksched.print_type_dict(typefile);

    for (uint_fast32_t iz = 0; iz < refZvec.size(); ++iz) {

      ctm_warning("iz: %" PRIuFAST32, iz);
      ctm_warning("nmax: %" PRIuFAST32, Tmatrices[iz]->get_nmax());

      assert_condition(Tmatrices[iz]->get_nmax() < maximum_order - 1);

      const float_type qext = Tmatrices[iz]->get_extinction_coefficient();
      const float_type qsca = Tmatrices[iz]->get_scattering_coefficient();
      const float_type walb = -qsca / qext;

      assert_values_equal_rel(double(qext), double(refQextvec[iz]), 1.e-5);
      assert_values_equal_rel(double(qsca), double(refQscavec[iz]), 1.e-5);
      assert_values_equal_rel(double(walb), double(refwalbvec[iz]), 1.e-5);

      // compare the result with the reference
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          assert_values_equal_rel(double((*Zmatrices[iz])(i, j)),
                                  double(refZvec[iz](i, j)), 0.1);
        }
      }
    }

    clear_vector(auxspace);
    clear_vector(quadrature_points);
    clear_vector(wignerdm0);
    clear_vector(converged_sizes);
    clear_vector(Tmatrices);
    clear_vector(Qtasks);
    clear_vector(Zmatrices);
    clear_vector(geometries);
    clear_vector(interactions);
    clear_vector(m0tasks);
    clear_vector(wignerdmalls);
    clear_vector(malltasks);
  }

  return 0;
}
