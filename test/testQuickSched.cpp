/**
 * @file testQuickSched.cpp
 *
 * @brief Unit test for the QuickSched library dependency.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Error.hpp"
#include "GaussBasedResources.hpp"
#include "NBasedResources.hpp"
#include "QuickSchedWrapper.hpp"
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
   * @brief MatrixMultiplicationTask
   *
   * @param i
   * @param j
   * @param m
   * @param k
   * @param a
   * @param b
   * @param c
   */
  MatrixMultiplicationTask(const uint_fast32_t i, const uint_fast32_t j,
                           const uint_fast32_t m, const uint_fast32_t k,
                           double *a, double *b, double *c)
      : _d{i, j}, _m(m), _k(k), _a(a), _b(b), _c(c) {}

  virtual ~MatrixMultiplicationTask() {}

  /**
   * @brief execute
   */
  virtual void execute() {
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
    QuickSched quicksched(4);

    NBasedResources nfactors(200);
    quicksched.register_resource(nfactors);
    quicksched.register_task(nfactors);
    quicksched.link_task_and_resource(nfactors, nfactors, true);

    std::vector<GaussBasedResources *> gaussfactors(100, nullptr);
    std::vector<WignerDResources *> wignerfactors(100, nullptr);
    for (uint_fast32_t ig = 0; ig < 100; ++ig) {
      gaussfactors[ig] = new GaussBasedResources(ig + 20);
      quicksched.register_resource(*gaussfactors[ig]);
      quicksched.register_task(*gaussfactors[ig]);
      quicksched.link_task_and_resource(*gaussfactors[ig], *gaussfactors[ig],
                                        true);

      wignerfactors[ig] = new WignerDResources(100, ig + 20, *gaussfactors[ig]);
      quicksched.register_resource(*wignerfactors[ig]);
      quicksched.register_task(*wignerfactors[ig]);
      quicksched.link_task_and_resource(*wignerfactors[ig], *wignerfactors[ig],
                                        true);
      quicksched.link_task_and_resource(*wignerfactors[ig], *gaussfactors[ig],
                                        false);
      quicksched.link_tasks(*gaussfactors[ig], *wignerfactors[ig]);
    }

    quicksched.execute_tasks(4);

    std::ofstream taskfile("test_quicksched_tasks.txt");
    taskfile << "# thread\tstart\tend\ttype\n";
    quicksched.print_task(nfactors, taskfile);
    for (uint_fast32_t ig = 0; ig < 100; ++ig) {
      quicksched.print_task(*gaussfactors[ig], taskfile);
      delete gaussfactors[ig];
      quicksched.print_task(*wignerfactors[ig], taskfile);
      delete wignerfactors[ig];
    }
    std::ofstream typefile("test_quicksched_types.txt");
    typefile << "# type\tlabel\n";
    quicksched.print_type_dict(typefile);
  }

  return 0;
}
