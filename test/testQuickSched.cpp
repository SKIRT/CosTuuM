/**
 * @file testQuickSched.cpp
 *
 * @brief Unit test for the QuickSched library dependency.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Error.hpp"
#include "GaussBasedResources.hpp"
#include "InteractionResource.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "TMatrixResource.hpp"
#include "Utilities.hpp"
#include "WignerDResources.hpp"

#include <cinttypes>
#include <fstream>
#include <vector>

/**
 * @brief Matrix multiplication kernel.
 *
 * @param m Number of rows in the first matrix.
 * @param n Number of columns in the first matrix and rows in the second matrix.
 * @param k Number of columns in the second matrix.
 * @param a First matrix.
 * @param lda Number of rows in storage layout of first matrix.
 * @param b Second matrix.
 * @param ldb Number of rows in storage layout of second matrix.
 * @param c Result matrix.
 * @param ldc Number of rows in storage layout of result matrix.
 */
void matmul(const uint_fast32_t m, const uint_fast32_t n, const uint_fast32_t k,
            const double *a, const uint_fast32_t lda, const double *b,
            const uint_fast32_t ldb, double *c, const uint_fast32_t ldc) {

  for (uint_fast32_t ii = 0; ii < m; ++ii)
    for (uint_fast32_t jj = 0; jj < n; ++jj) {
      double acc = 0.0;
      for (uint_fast32_t kk = 0; kk < k; ++kk)
        acc += a[ii + lda * kk] * b[kk + ldb * jj];
      c[ii + ldc * jj] += acc;
    }
}

/**
 * @brief Task that multiplies a single block in the result matrix.
 */
class MatrixMultiplicationTask : public Task {
private:
  /*! @brief Offsets of a block to calculate in the result matrix. */
  const uint_fast32_t _d[2];

  /*! @brief Number of rows in the first matrix. */
  const uint_fast32_t _m;

  /*! @brief Number of columns in the second matrix. */
  const uint_fast32_t _k;

  /*! @brief First matrix. */
  const double *_a;

  /*! @brief Second matrix. */
  const double *_b;

  /*! @brief Result matrix. */
  double *_c;

public:
  /**
   * @brief Constructor.
   *
   * @param i Row index.
   * @param j Column index.
   * @param m Number of rows in the first matrix.
   * @param k Number of columns in the second matrix.
   * @param a Matrix A.
   * @param b Matrix B.
   * @param c Result matrix.
   */
  MatrixMultiplicationTask(const uint_fast32_t i, const uint_fast32_t j,
                           const uint_fast32_t m, const uint_fast32_t k,
                           double *a, double *b, double *c)
      : _d{i, j}, _m(m), _k(k), _a(a), _b(b), _c(c) {}

  virtual ~MatrixMultiplicationTask() {}

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {
    matmul(32, 32, _k * 32, &_a[_d[0] * 32], _m * 32, &_b[_k * 32 * _d[1] * 32],
           _k * 32, &_c[_d[0] * 32 + _m * 32 * _d[1] * 32], _m * 32);
  }
};

/**
 * @brief Unit test for the QuickSched library dependency.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// basic test (copied from QuickSched itself)
  {
    double *a, *b, *c, irm = 1.0 / RAND_MAX;
    const uint_fast32_t m = 30;
    const uint_fast32_t n = 30;
    const uint_fast32_t k = 30;

    a = new double[32 * 32 * m * k];
    b = new double[32 * 32 * k * n];
    c = new double[32 * 32 * m * n];

    /* Fill the matrices. */
    for (uint_fast32_t i = 0; i < m * k * 32 * 32; ++i) {
      a[i] = rand() * irm;
    }
    for (uint_fast32_t i = 0; i < k * n * 32 * 32; ++i) {
      b[i] = rand() * irm;
    }
    for (uint_fast32_t i = 0; i < k * n * 32 * 32; ++i) {
      c[i] = 0.;
    }

    std::vector<MatrixMultiplicationTask> tasks;
    // the resources are simple dummies in this case
    std::vector<Resource> resources(m * n);
    for (uint_fast32_t i = 0; i < m; ++i) {
      for (uint_fast32_t j = 0; j < n; ++j) {
        tasks.push_back(MatrixMultiplicationTask(i, j, m, k, a, b, c));
      }
    }

    QuickSched quicksched(4);
    /* Build a task for each tile of the matrix c. */
    for (uint_fast32_t i = 0; i < tasks.size(); ++i) {
      quicksched.register_resource(resources[i]);
      quicksched.register_task(tasks[i]);
      quicksched.link_task_and_resource(tasks[i], resources[i], true);
    }

    quicksched.execute_tasks(4);

    delete[] a;
    delete[] b;
    delete[] c;
  }

  {
    const uint_fast32_t auxsize = 100;
    const uint_fast32_t maxgauss = 200;
    const uint_fast32_t maximum_order = 100;

    size_t size = 0;

    const size_t nbasedresourcessize =
        NBasedResources::get_memory_size(maximum_order);
    ctm_warning("NBasedResources: %s",
                Utilities::human_readable_bytes(nbasedresourcessize).c_str());
    size += nbasedresourcessize;

    const size_t auxiliarysize =
        TMatrixAuxiliarySpace::get_memory_size(maximum_order);
    ctm_warning("TMatrixAuxiliarySpace: %s",
                Utilities::human_readable_bytes(auxiliarysize).c_str());
    size += auxsize * auxiliarysize;

    const size_t gaussresourcesize =
        GaussBasedResources::get_memory_size(maxgauss);
    ctm_warning("GaussBasedResources: %s",
                Utilities::human_readable_bytes(gaussresourcesize).c_str());
    size += maxgauss * gaussresourcesize;

    const size_t wignerm0resourcesize =
        WignerDm0Resources::get_memory_size(maximum_order, maxgauss);
    ctm_warning("WignerDResources: %s",
                Utilities::human_readable_bytes(wignerm0resourcesize).c_str());
    size += maxgauss * wignerm0resourcesize;

    const size_t geometryresourcesize =
        ParticleGeometryResource::get_memory_size(maxgauss);
    ctm_warning("ParticleGeometryResource: %s",
                Utilities::human_readable_bytes(geometryresourcesize).c_str());
    size += maxgauss * geometryresourcesize;

    const size_t interactionresourcesize =
        InteractionResource::get_memory_size(maximum_order, maxgauss);
    ctm_warning(
        "InteractionResource: %s",
        Utilities::human_readable_bytes(interactionresourcesize).c_str());
    size += maxgauss * interactionresourcesize;

    const size_t tmatrixresourcesize =
        TMatrixResource::get_memory_size(maximum_order);
    ctm_warning("TMatrixResource: %s",
                Utilities::human_readable_bytes(tmatrixresourcesize).c_str());
    size += maxgauss * tmatrixresourcesize;

    ctm_warning("Total: %s", Utilities::human_readable_bytes(size).c_str());

    QuickSched quicksched(4);

    NBasedResources nfactors(maximum_order);
    quicksched.register_resource(nfactors);
    quicksched.register_task(nfactors);
    nfactors.link_resources(quicksched);

    ConvergedSizeResources converged_size;
    quicksched.register_resource(converged_size);

    TMatrixAuxiliarySpaceManager aux_manager(4, maximum_order);

    std::vector<GaussBasedResources *> gaussfactors(maxgauss, nullptr);
    std::vector<WignerDm0Resources *> wignerm0factors(maxgauss, nullptr);
    std::vector<ParticleGeometryResource *> geometryfactors(maxgauss, nullptr);
    std::vector<InteractionResource *> interactionfactors(maxgauss, nullptr);
    std::vector<InteractionTask *> interaction_tasks(maxgauss, nullptr);
    std::vector<TMatrixResource *> tmatrices(maxgauss, nullptr);
    std::vector<TMatrixM0Task *> tmatrixm0tasks(maxgauss, nullptr);
    for (uint_fast32_t ig = 0; ig < maxgauss; ++ig) {
      gaussfactors[ig] = new GaussBasedResources(ig + 20);
      quicksched.register_resource(*gaussfactors[ig]);
      quicksched.register_task(*gaussfactors[ig]);
      gaussfactors[ig]->link_resources(quicksched);

      wignerm0factors[ig] =
          new WignerDm0Resources(maximum_order, ig + 20, *gaussfactors[ig]);
      quicksched.register_resource(*wignerm0factors[ig]);
      quicksched.register_task(*wignerm0factors[ig]);
      wignerm0factors[ig]->link_resources(quicksched);
      quicksched.link_tasks(*gaussfactors[ig], *wignerm0factors[ig]);

      geometryfactors[ig] =
          new ParticleGeometryResource(10., 0.5, ig + 20, *gaussfactors[ig]);
      quicksched.register_resource(*geometryfactors[ig]);
      quicksched.register_task(*geometryfactors[ig]);
      geometryfactors[ig]->link_resources(quicksched);
      quicksched.link_tasks(*gaussfactors[ig], *geometryfactors[ig]);

      interactionfactors[ig] = new InteractionResource(
          100., std::complex<float_type>(1., 0.02), maximum_order, ig + 20);
      quicksched.register_resource(*interactionfactors[ig]);
      interaction_tasks[ig] =
          new InteractionTask(maximum_order, ig + 20, *geometryfactors[ig],
                              *interactionfactors[ig]);
      quicksched.register_task(*interaction_tasks[ig]);
      interaction_tasks[ig]->link_resources(quicksched);
      quicksched.link_tasks(*geometryfactors[ig], *interaction_tasks[ig]);

      tmatrices[ig] = new TMatrixResource(maximum_order);
      quicksched.register_resource(tmatrices[ig]->get_m_resource(0));

      tmatrixm0tasks[ig] = new TMatrixM0Task(
          1.e-4, 50, ig + 20, nfactors, *gaussfactors[ig], *geometryfactors[ig],
          *interactionfactors[ig], *wignerm0factors[ig], aux_manager,
          *tmatrices[ig], converged_size, tmatrices[ig]->get_m_resource(0));
      quicksched.register_task(*tmatrixm0tasks[ig]);
      tmatrixm0tasks[ig]->link_resources(quicksched);
      quicksched.link_tasks(nfactors, *tmatrixm0tasks[ig]);
      quicksched.link_tasks(*interaction_tasks[ig], *tmatrixm0tasks[ig]);
      quicksched.link_tasks(*wignerm0factors[ig], *tmatrixm0tasks[ig]);
    }

    quicksched.execute_tasks(4);

    std::ofstream taskfile("test_quicksched_tasks.txt");
    taskfile << "# thread\tstart\tend\ttype\n";
    quicksched.print_task(nfactors, taskfile);
    for (uint_fast32_t ig = 0; ig < maxgauss; ++ig) {
      quicksched.print_task(*gaussfactors[ig], taskfile);
      delete gaussfactors[ig];
      quicksched.print_task(*wignerm0factors[ig], taskfile);
      delete wignerm0factors[ig];
      quicksched.print_task(*geometryfactors[ig], taskfile);
      delete geometryfactors[ig];
      quicksched.print_task(*interaction_tasks[ig], taskfile);
      delete interactionfactors[ig];
      delete interaction_tasks[ig];
      quicksched.print_task(*tmatrixm0tasks[ig], taskfile);
      delete tmatrices[ig];
      delete tmatrixm0tasks[ig];
    }
    std::ofstream typefile("test_quicksched_types.txt");
    typefile << "# type\tlabel\n";
    quicksched.print_type_dict(typefile);
  }

  return 0;
}
