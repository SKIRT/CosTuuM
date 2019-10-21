/**
 * @file testQuickSched.cpp
 *
 * @brief Unit test for the QuickSched library dependency.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Error.hpp"
#include "quicksched.h"

#include <cinttypes>

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
class MatrixMultiplicationTask {
private:
  /*! @brief Offsets of a block to calculate in the result matrix. */
  const uint_fast32_t _d[2];

  /*! @brief Number of rows in the first matrix. */
  const uint_fast32_t _m;

  /*! @brief Number of columns in the first matrix and rows in the second
   *  matrix. */
  const uint_fast32_t _n;

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
   * @param n
   * @param k
   * @param a
   * @param b
   * @param c
   */
  MatrixMultiplicationTask(const uint_fast32_t i, const uint_fast32_t j,
                           const uint_fast32_t m, const uint_fast32_t n,
                           const uint_fast32_t k, double *a, double *b,
                           double *c)
      : _d{i, j}, _m(m), _n(n), _k(k), _a(a), _b(b), _c(c) {}

  /**
   * @brief execute
   */
  void execute() {
    matmul(32, 32, _k * 32, &_a[_d[0] * 32], _m * 32, &_b[_k * 32 * _d[1] * 32],
           _k * 32, &_c[_d[0] * 32 + _m * 32 * _d[1] * 32], _m * 32);
  }
};

/**
 * @brief Execute the given void-wrapped task.
 *
 * @param data Void wrapped task data.
 * @tparam TASK_TYPE Type of task to execute.
 */
template <typename TASK_TYPE> inline void execute(void *data) {
  TASK_TYPE *task = static_cast<TASK_TYPE *>(data);
  task->execute();
}

/**
 * @brief Runner function to pass on to the scheduler.
 *
 * @param type Type of task being run.
 * @param data Additional data for the task.
 */
void runner(int type, void *data) {

  switch (type) {
  case 1:
    execute<MatrixMultiplicationTask>(data);
    break;
  default:
    ctm_error("Unknown task type: %i", type);
  }
}

/**
 * @brief Unit test for the QuickSched library dependency.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const int_fast32_t nr_threads = 4;

  qsched_task_t tid;
  qsched_res_t rid;
  struct qsched s;
  double *a, *b, *c, irm = 1.0 / RAND_MAX;
  const uint_fast32_t m = 4;
  const uint_fast32_t n = 4;
  const uint_fast32_t k = 4;

  /* Tell the user something about the test. */
  message("computing a tiled matrix multiplication of the form "
          "C_ij = A_i: * B_:j, with a single task per C_ij.");

  /* Init the sched. */
  bzero(&s, sizeof(struct qsched));
  qsched_init(&s, nr_threads, qsched_flag_none);

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

  /* Build a task for each tile of the matrix c. */
  for (uint_fast32_t i = 0; i < m; ++i)
    for (uint_fast32_t j = 0; j < n; ++j) {
      MatrixMultiplicationTask task(i, j, m, n, k, a, b, c);
      rid = qsched_addres(&s, qsched_owner_none, qsched_res_none);
      tid = qsched_addtask(&s, 1, task_flag_none, &task,
                           sizeof(MatrixMultiplicationTask), 1);
      qsched_addlock(&s, tid, rid);
    }

  qsched_run(&s, nr_threads, runner);

  /* Clean up. */
  qsched_free(&s);

  delete[] a;
  delete[] b;
  delete[] c;

  return 0;
}
