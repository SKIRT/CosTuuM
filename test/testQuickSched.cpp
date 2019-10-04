/**
 * @file testQuickSched.cpp
 *
 * @brief Unit test for the QuickSched library dependency.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

extern "C" {
#include "quicksched.h"
}

/**
 * @brief Matrix multiplication kernel.
 */

void matmul(int m, int n, int k, const double *a, int lda, const double *b,
            int ldb, double *c, int ldc) {

  int ii, jj, kk;
  double acc;

  for (ii = 0; ii < m; ii++)
    for (jj = 0; jj < n; jj++) {
      for (acc = 0.0, kk = 0; kk < k; kk++)
        acc += a[ii + lda * kk] * b[kk + ldb * jj];
      c[ii + ldc * jj] += acc;
    }
}

struct task_data {
  int d[2];
  int m;
  int n;
  int k;
  double *a;
  double *b;
  double *c;
};

/* Runner function to pass to the scheduler. */
void runner(int type, void *data) {

  /* Decode the task data. */
  struct task_data *tdata = (struct task_data *)data;
  const int k = tdata->k;
  const int m = tdata->m;
  const double *a = tdata->a;
  const double *b = tdata->b;
  double *c = tdata->b;
  const int *d = tdata->d;

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

  const int nr_threads = 4;

  int i, j;
  qsched_task_t tid;
  qsched_res_t rid;
  struct qsched s;
  double *a, *b, *c, irm = 1.0 / RAND_MAX;
  int m = 4;
  int n = 4;
  int k = 4;

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
  for (i = 0; i < m * k * 32 * 32; i++) {
    a[i] = rand() * irm;
  }
  for (i = 0; i < k * n * 32 * 32; i++) {
    b[i] = rand() * irm;
  }
  for (i = 0; i < k * n * 32 * 32; i++) {
    c[i] = 0.;
  }

  /* Build a task for each tile of the matrix c. */
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
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
