/**
 * @file testQuickSched.cpp
 *
 * @brief Unit test for the QuickSched library dependency.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "quicksched.h"

/**
 * @brief Matrix multiplication kernel.
 */

void matmul(int m, int n, int k, double *a, int lda, double *b, int ldb,
            double *c, int ldc) {

  int ii, jj, kk;
  double acc;

  // message( "matmul: m=%i, n=%i, k=%i, lda=%i, ldb=%i, ldc=%i." ,
  //     m , n , k , lda , ldb , ldc ); fflush(stdout);

  for (ii = 0; ii < m; ii++)
    for (jj = 0; jj < n; jj++) {
      for (acc = 0.0, kk = 0; kk < k; kk++)
        acc += a[ii + lda * kk] * b[kk + ldb * jj];
      c[ii + ldc * jj] += acc;
    }
}

/* Runner function to pass to the scheduler. */
void runner(int type, void *data) {

  /* Decode the task data. */
  int *d = (int *)data;

  /* Decode and execute the task. */
  switch (type) {
  case 1:
    matmul(32, 32, k * 32, &a[d[0] * 32], m * 32, &b[k * 32 * d[1] * 32],
           k * 32, &c[d[0] * 32 + m * 32 * d[1] * 32], m * 32);
    break;
  default:
    error("Unknown task type.");
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

  int i, j, data[2];
  qsched_task_t tid;
  qsched_res_t rid;
  struct qsched s;
  double *a, *b, *c, *res, irm = 1.0 / RAND_MAX;

  /* Tell the user something about the test. */
  message("computing a tiled matrix multiplication of the form "
          "C_ij = A_i: * B_:j, with a single task per C_ij.");

  /* Init the sched. */
  bzero(&s, sizeof(struct qsched));
  qsched_init(&s, nr_threads, qsched_flag_none);

  /* Allocate the matrices. */
  if ((a = (double *)malloc(sizeof(double) * m * k * 32 * 32)) == NULL ||
      (b = (double *)malloc(sizeof(double) * k * n * 32 * 32)) == NULL ||
      (c = (double *)malloc(sizeof(double) * m * n * 32 * 32)) == NULL ||
      (res = (double *)malloc(sizeof(double) * m * n * 32 * 32)) == NULL)
    error("Failed to allocate matrices.");

  /* Fill the matrices. */
  for (i = 0; i < m * k * 32 * 32; i++)
    a[i] = rand() * irm;
  for (i = 0; i < k * n * 32 * 32; i++)
    b[i] = rand() * irm;
  bzero(c, sizeof(double) * m * n * 32 * 32);
  bzero(res, sizeof(double) * m * n * 32 * 32);

  /* Build a task for each tile of the matrix c. */
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      data[0] = i;
      data[1] = j;
      rid = qsched_addres(&s, qsched_owner_none, qsched_res_none);
      tid = qsched_addtask(&s, 1, task_flag_none, data, 2 * sizeof(int), 1);
      qsched_addlock(&s, tid, rid);
    }

  qsched_run(&s, nr_threads, runner);

  /* Clean up. */
  qsched_free(&s);

  return 0;
}
