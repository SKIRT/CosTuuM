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
 * @brief Data struct used to store data for the test tasks.
 */
struct task_data {
  /*! @brief Offsets of a block to calculate in the result matrix. */
  uint_fast32_t d[2];
  /*! @brief Number of rows in the first matrix. */
  uint_fast32_t m;
  /*! @brief Number of columns in the first matrix and rows in the second
   *  matrix. */
  uint_fast32_t n;
  /*! @brief Number of columns in the second matrix. */
  uint_fast32_t k;
  /*! @brief First matrix. */
  double *a;
  /*! @brief Second matrix. */
  double *b;
  /*! @brief Result matrix. */
  double *c;
};

/**
 * @brief Runner function to pass on to the scheduler.
 *
 * @param type Type of task being run.
 * @param data Additional data for the task.
 */
void runner(int type, void *data) {

  /* Decode the task data. */
  struct task_data *tdata = (struct task_data *)data;
  const uint_fast32_t k = tdata->k;
  const uint_fast32_t m = tdata->m;
  const double *a = tdata->a;
  const double *b = tdata->b;
  double *c = tdata->b;
  const uint_fast32_t *d = tdata->d;

  /* Decode and execute the task. */
  switch (type) {
  case 1:
    matmul(32, 32, k * 32, &a[d[0] * 32], m * 32, &b[k * 32 * d[1] * 32],
           k * 32, &c[d[0] * 32 + m * 32 * d[1] * 32], m * 32);
    break;
  default:
    ctm_error("Unknown task type.");
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
      struct task_data tdata;
      tdata.d[0] = i;
      tdata.d[1] = j;
      tdata.m = m;
      tdata.n = n;
      tdata.k = k;
      tdata.a = a;
      tdata.b = b;
      tdata.c = c;
      rid = qsched_addres(&s, qsched_owner_none, qsched_res_none);
      tid = qsched_addtask(&s, 1, task_flag_none, &tdata,
                           sizeof(struct task_data), 1);
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
