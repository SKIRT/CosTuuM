/*******************************************************************************
 * This file is part of QuickSched.
 * Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "config.h"

/* Standard includes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* OpenMP headers, only if available. */
#ifdef QUICKSCHED_HAVE_OPENMP
#include <omp.h>
#endif

// /* Pthread headers, only if available. */
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

/* Local includes. */
#include "atomic.h"
#include "cycle.h"
#include "error.h"
#include "lock.h"
#include "qsched.h"
#include "queue.h"
#include "res.h"
#include "task.h"

/** Timer names. */
char *qsched_timer_names[qsched_timer_count] = {"queue",   "qlock", "lock",
                                                "gettask", "done",  "prepare"};

/**
 * @brief Change the owner of a resource.
 *
 * @param s Pointer to the #qsched.
 * @param res Resource handle.
 * @param owner The ID of the new owner.
 */

void qsched_res_own(struct qsched *s, qsched_res_t res, int owner) {

  /* Force the new ownership. */
  s->res[res].owner = owner;
}

/**
 * @brief Make sure that the scheduler has enough memory
 *        allocated for the tasks, resources, dependencies,
 *        locks, uses, and data.
 *
 * @param s Pointer to the #qsched.
 * @param nr_tasks Minimum number of tasks.
 * @param nr_res Minimum number of reources.
 * @param nr_deps Minimum number of dependencies.
 * @param nr_locks Minimum number of locks.
 * @param nr_uses Minimum number of uses.
 * @param size_data Mnimum size of task data.
 *
 * This procedure pre-allocates resources in the sched. Although
 * These buffers are grown automatically, this is not possible
 * when dynamic tasks are spawned (see #qsched_addtask_dynamic),
 * and insufficient buffers my cause dynamic task allocation to fail!
 *
 * Keep in mind that the individual blocks of task data are
 * rounded up to 16 bytes for alignment.
 * The number of dependencies, locks, and uses refers to the
 * number of times #qsched_addunlock, #qsched_addlock, and
 * #qsched_adduse are called, respectively.
 */

void qsched_ensure(struct qsched *s, int nr_tasks, int nr_res, int nr_deps,
                   int nr_locks, int nr_uses, int size_data) {

  int dirty = 0;

  /* Re-allocate tasks? */
  if (s->size < nr_tasks) {
    dirty = 1;
    struct task *tasks_new;
    if ((tasks_new = (struct task *)malloc(sizeof(struct task) * nr_tasks)) ==
        NULL)
      quicksched_error("Failed to allocate new task buffer.");
    memcpy(tasks_new, s->tasks, sizeof(struct task) * s->count);
    free(s->tasks);
    s->tasks = tasks_new;
    s->size = nr_tasks;
  }

  /* Re-allocate resources? */
  if (s->size_res < nr_res) {
    dirty = 1;
    struct res *res_new;
    if ((res_new = (struct res *)malloc(sizeof(struct res) * nr_res)) == NULL)
      quicksched_error("Failed to allocate new res buffer.");
    memcpy(res_new, s->res, sizeof(struct res) * s->count_res);
    free(s->res);
    s->res = res_new;
    s->size_res = nr_res;
  }

  /* Re-allocate dependencies? */
  if (s->size_deps < nr_deps) {
    dirty = 1;
    qsched_task_t *deps_new, *deps_key_new;
    if ((deps_new = (qsched_task_t *)malloc(sizeof(qsched_task_t) * nr_deps)) ==
            NULL ||
        (deps_key_new =
             (qsched_task_t *)malloc(sizeof(qsched_task_t) * nr_deps)) == NULL)
      quicksched_error("Failed to allocate new deps buffer.");
    memcpy(deps_new, s->deps, sizeof(qsched_task_t) * s->count_deps);
    memcpy(deps_key_new, s->deps_key, sizeof(qsched_task_t) * s->count_deps);
    free(s->deps);
    free(s->deps_key);
    s->deps = deps_new;
    s->deps_key = deps_key_new;
    s->size_deps = nr_deps;
  }

  /* Re-allocate locks? */
  if (s->size_locks < nr_locks) {
    dirty = 1;
    qsched_res_t *locks_new;
    qsched_task_t *locks_key_new;
    if ((locks_new = (qsched_res_t *)malloc(sizeof(qsched_res_t) * nr_locks)) ==
            NULL ||
        (locks_key_new =
             (qsched_task_t *)malloc(sizeof(qsched_task_t) * nr_locks)) == NULL)
      quicksched_error("Failed to allocate new locks buffer.");
    memcpy(locks_new, s->locks, sizeof(qsched_res_t) * s->count_locks);
    memcpy(locks_key_new, s->locks_key, sizeof(qsched_task_t) * s->count_locks);
    free(s->locks);
    free(s->locks_key);
    s->locks = locks_new;
    s->locks_key = locks_key_new;
    s->size_locks = nr_locks;
  }

  /* Re-allocate uses? */
  if (s->size_uses < nr_uses) {
    dirty = 1;
    qsched_res_t *uses_new;
    qsched_task_t *uses_key_new;
    if ((uses_new = (qsched_res_t *)malloc(sizeof(qsched_res_t) * nr_uses)) ==
            NULL ||
        (uses_key_new =
             (qsched_task_t *)malloc(sizeof(qsched_task_t) * nr_uses)) == NULL)
      quicksched_error("Failed to allocate new uses buffer.");
    memcpy(uses_new, s->uses, sizeof(qsched_res_t) * s->count_uses);
    memcpy(uses_key_new, s->uses_key, sizeof(qsched_task_t) * s->count_uses);
    free(s->uses);
    free(s->uses_key);
    s->uses = uses_new;
    s->uses_key = uses_key_new;
    s->size_uses = nr_uses;
  }

  /* Re-allocate resources? */
  if (s->size_data < size_data) {
    dirty = 1;
    char *data_new;
    if ((data_new = (char *)malloc(size_data)) == NULL)
      quicksched_error("Failed to allocate new data buffer.");
    memcpy(data_new, s->data, s->count_data);
    free(s->data);
    s->data = data_new;
    s->size_data = size_data;
  }

  /* Mark scheduler if dirty. */
  if (dirty)
    s->flags |= qsched_flag_dirty;
}

/**
 * @brief Add a task to the scheduler on the fly.
 *
 * @param s Pointer to the #qsched.
 * @param type Task type.
 * @param flags Task flags.
 * @param data Task data.
 * @param data_size Size, in bytes, of task data.
 * @param cost Relative task computational cost.
 * @param locks Array of #qsched_res_t locked by this task.
 * @param nr_locks Number of locked resources.
 * @param uses Array of #qsched_res_t used by this task.
 * @param nr_uses Number of used resources.
 */

void qsched_addtask_dynamic(struct qsched *s, int type, unsigned int flags,
                            void *data, int data_size, int cost,
                            qsched_res_t *locks, int nr_locks,
                            qsched_res_t *uses, int nr_uses) {

  int k, tid, ind, data_size2;
  struct task *t;

  /* Allocate a new task. */
  if ((tid = atomic_inc(&s->count)) >= s->size)
    quicksched_error("Task buffer overflow.");
  t = &s->tasks[tid];

  /* Set the task data. */
  t->type = type;
  t->flags = flags;
  t->weight = t->cost = cost;
  t->nr_unlocks = 0;
  t->unlocks = NULL;

  /* Add the data. */
  data_size2 = (data_size + (qsched_data_round - 1)) & ~(qsched_data_round - 1);
  if ((ind = atomic_add(&s->count_data, data_size2)) + data_size2 >
      s->size_data)
    quicksched_error("Data buffer overflow.");
  memcpy(&s->data[ind], data, data_size);
  t->data = ind;

  /* Add the locks. */
  if ((ind = atomic_add(&s->count_locks, nr_locks)) + nr_locks > s->size_locks)
    quicksched_error("Locks buffer overflow.");
  memcpy(&s->locks[ind], locks, sizeof(qsched_res_t) * nr_locks);
  for (k = 0; k < nr_locks; k++)
    s->locks_key[ind + k] = tid;
  t->locks = &s->locks[ind];

  /* Add the uses. */
  if ((ind = atomic_add(&s->count_uses, nr_uses)) + nr_uses > s->size_uses)
    quicksched_error("uses buffer overflow.");
  memcpy(&s->uses[ind], uses, sizeof(qsched_res_t) * nr_uses);
  for (k = 0; k < nr_uses; k++)
    s->uses_key[ind + k] = tid;
  t->uses = &s->uses[ind];

  /* The task is now ready to run, submit it. */
  t->wait = 0;
  qsched_enqueue(s, t);
}

/**
 * @brief Clear the tasks and resources in a scheduler.
 *
 * @param s Pointer to the #qsched.
 */

void qsched_reset(struct qsched *s) {

  /* Simply clear the counts, leave the buffers allocated. */
  s->count = 0;
  s->waiting = 0;
  s->count_data = 0;
  s->count_deps = 0;
  s->count_locks = 0;
  s->count_uses = 0;
  s->count_res = 0;

/* Clear the timers. */
#ifdef TIMERS
  bzero(s->timers, sizeof(ticks) * qsched_timer_count);
#endif
}

/**
 * @brief Execute all the tasks in the current scheduler using
 *        OpenMP.
 *
 * @param s Pointer to the #qsched.
 * @param nr_threads Number of threads to use.
 * @param fun User-supplied function that will be called with the
 *        task type and a pointer to the task data.
 *
 * This function is only available if QuickSched was compiled with
 * OpenMP support.
 */

void qsched_run_openmp(struct qsched *s, int nr_threads, qsched_funtype fun) {

#if defined(QUICKSCHED_HAVE_OPENMP)

  /* Prepare the scheduler. */
  qsched_prepare(s);

/* Parallel loop. */
#pragma omp parallel num_threads(nr_threads)
  {
    /* Local variable. */
    struct task *t;

    /* Get the ID of the current thread. */
    int qid = omp_get_thread_num() % s->nr_queues;

    /* Loop as long as there are tasks. */
    while ((t = qsched_gettask(s, qid)) != NULL) {

      /* Call the user-supplied function on the task with its data. */
      fun(t->type, &s->data[t->data]);

      /* Mark that task as done. */
      qsched_done(s, t);

    } /* loop as long as there are tasks. */

  } /* parallel loop. */

#else
  error("QuickSched was not compiled with OpenMP support.");
#endif
}

/**
 * @brief Barrier function when running with pthreads.
 *
 * @param s Pointer to the #qsched.
 * @param tid ID of the calling thread. This is needed in case
 *        we don't want to launch all threads.
 */

void qsched_barrier_wait(struct qsched *s, int tid) {

#if defined(HAVE_PTHREAD)
  /* First, get the barrier mutex. */
  if (pthread_mutex_lock(&s->barrier_mutex) != 0)
    quicksched_error("Failed to get barrier mutex.");

  /* The callee's thread stops running. */
  s->barrier_running -= 1;

  /* If all threads are in, send a signal... */
  if (s->barrier_running == 0)
    if (pthread_cond_broadcast(&s->barrier_cond) != 0)
      quicksched_error("Failed to broadcast barrier full condition.");

  /* Wait for the barrier to open. */
  while (s->barrier_count == 0 || tid >= s->barrier_launchcount)
    if (pthread_cond_wait(&s->barrier_cond, &s->barrier_mutex) != 0)
      quicksched_error("Eror waiting for barrier to close.");

  /* This thread is leaving, decrease the barrier count, increase
     the number of threads running. */
  s->barrier_count -= 1;
  s->barrier_running += 1;

  /* If I'm the last one out, signal the condition again. */
  if (s->barrier_count == 0)
    if (pthread_cond_broadcast(&s->barrier_cond) != 0)
      quicksched_error("Failed to broadcast empty barrier condition.");

  /* Last but not least, release the mutex. */
  if (pthread_mutex_unlock(&s->barrier_mutex) != 0)
    quicksched_error("Failed to get unlock the barrier mutex.");
#endif
}

/**
 * @brief Launch a given number of threads and wait for them to finish.
 *
 * @param s Pointer to the #qsched.
 * @param nr_threads Number of threads to let through the barrier.
 *
 * Lets the specified number of threads through the barrier. Note that
 * this assumes that at least that many threads are waiting there,
 * otherwise this function will hang indefinitely.
 */
void qsched_launch_threads(struct qsched *s, int nr_threads) {

#if defined(HAVE_PTHREAD)
  /* Wait for all the runners to have entered the barrier. */
  while (s->barrier_running)
    if (pthread_cond_wait(&s->barrier_cond, &s->barrier_mutex) != 0)
      quicksched_error("Error while waiting for barrier.");

  /* Cry havoc and let loose the dogs of war. */
  s->barrier_count = nr_threads;
  s->barrier_launchcount = nr_threads;
  if (pthread_cond_broadcast(&s->barrier_cond) != 0)
    quicksched_error("Failed to broadcast barrier open condition.");

  /* Lean back and wait for the runners to come home. */
  while (s->barrier_count || s->barrier_running)
    if (pthread_cond_wait(&s->barrier_cond, &s->barrier_mutex) != 0)
      quicksched_error("Error while waiting for barrier.");
#endif
}

void *qsched_pthread_run(void *in) {

  struct qsched_pthread_runner *r = (struct qsched_pthread_runner *)in;
  struct qsched *s = r->s;
  int tid = r->tid;
  struct task *t;

  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    qsched_barrier_wait(s, tid);

    /* If there is no function to execute, then just quit. */
    if (s->fun == NULL)
      pthread_exit(NULL);

    /* Loop as long as there are tasks. */
    while ((t = qsched_gettask(s, tid)) != NULL) {

      /* Call the user-supplied function on the task with its data. */
      s->fun(t->type, &s->data[t->data]);

      /* Mark that task as done. */
      qsched_done(s, t);

    } /* loop as long as there are tasks. */

  } /* main loop. */
}

/**
 * @brief Execute all the tasks in the current scheduler using
 *        pthreads.
 *
 * @param s Pointer to the #qsched.
 * @param nr_threads Number of threads to use.
 * @param fun User-supplied function that will be called with the
 *        task type and a pointer to the task data.
 *
 * This function is only available if QuickSched was compiled with
 * pthread support.
 */

void qsched_run_pthread(struct qsched *s, int nr_threads, qsched_funtype fun) {

#if defined(HAVE_PTHREAD)

  /* Prepare the scheduler. */
  qsched_prepare(s);

  /* Make sure we have enough threads. */
  if (nr_threads > s->runners_count) {

    /* Reallocate the threads array? */
    if (nr_threads > s->runners_size) {
      int runners_size_new = s->runners_size;
      while (runners_size_new < nr_threads)
        runners_size_new *= qsched_stretch;
      struct qsched_pthread_runner *runners_new;
      if ((runners_new = malloc(sizeof(struct qsched_pthread_runner) *
                                runners_size_new)) == NULL)
        quicksched_error("Failed to allocate new thread array.");
      memcpy(runners_new, s->runners, sizeof(pthread_t) * s->runners_count);
      free(s->runners);
      s->runners = runners_new;
      s->runners_size = runners_size_new;
    }

    /* Launch the missing threads. */
    for (int tid = s->runners_count; tid < nr_threads; tid++) {
      s->barrier_running += 1;
      s->runners[tid].tid = tid;
      s->runners[tid].s = s;
      if (pthread_create(&s->runners[tid].thread, NULL, qsched_pthread_run,
                         (void *)&s->runners[tid]) != 0)
        quicksched_error("Failed to create pthread.");
      s->runners_count += 1;
    }
  }

  /* Set the runner function. */
  s->fun = fun;

  /* Launch the threads. */
  qsched_launch_threads(s, nr_threads);

#else
  error("QuickSched was not compiled with pthread support.");
#endif
}

/**
 * @brief Execute all the tasks in the current scheduler using
 *        pthreads.
 *
 * @param s Pointer to the #qsched.
 * @param nr_threads Number of threads to use.
 * @param fun User-supplied function that will be called with the
 *        task type and a pointer to the task data.
 *
 * Depending on what type of threads QuickSched was compiled with,
 * i.e. OpenMP and/or pthreads, and if the #qsched_flag_spin or
 * #qsched_flag_pthread flags were set, this function calls
 * #qsched_run_openmp or #qsched_run_pthread respectively.
 */

void qsched_run(struct qsched *s, int nr_threads, qsched_funtype fun) {

  /* Force use of pthreads? */
  if (!QUICKSCHED_HAVE_OPENMP ||
      s->flags & (qsched_flag_yield | qsched_flag_pthread))
    qsched_run_pthread(s, nr_threads, fun);

  /* Otherwise, default to OpenMP. */
  else if (QUICKSCHED_HAVE_OPENMP)
    qsched_run_openmp(s, nr_threads, fun);

  else
    quicksched_error("QuickSched was not compiled with OpenMP or pthreads.");
}

/**
 * @brief Fetch the data pointer of a task.
 *
 * @param s Pointer to the #qsched.
 * @param t Pointer to the #task.
 */

void *qsched_getdata(struct qsched *s, struct task *t) {

  return &s->data[t->data];
}

/**
 * @brief Put the given task in the best possible queue.
 *
 * @param s Pointer to the #qsched.
 * @param t Pointer to the #task.
 */

void qsched_enqueue(struct qsched *s, struct task *t) {

  int j, qid, scores[s->nr_queues], oid;

  /* If this is a virtual task, just do its unlocks and leave. */
  if (t->flags & task_flag_virtual) {

    /* This task is done before it started. */
    t->tic = getticks();
    qsched_done(s, t);

  }

  /* Otherwise, find a home (queue) for it. */
  else {

    /* Init the scores for each queue. */
    for (j = 0; j < s->nr_queues; j++)
      scores[j] = 0;

    /* Loop over the locks and uses, and get their owners. */
    for (j = 0; j < t->nr_locks; j++)
      if ((oid = s->res[t->locks[j]].owner) != qsched_owner_none)
        scores[oid] += 1;
    for (j = 0; j < t->nr_uses; j++)
      if ((oid = s->res[t->uses[j]].owner) != qsched_owner_none)
        scores[oid] += 1;

    /* Find the queue with the highest score. */
    qid = 0;
    for (j = 1; j < s->nr_queues; j++)
      if (scores[j] > scores[qid] ||
          (scores[j] == scores[qid] &&
           s->queues[j].count < s->queues[qid].count))
        qid = j;

    /* Put the unlocked task in that queue. */
    queue_put(&s->queues[qid], s, t - s->tasks);
  }
}

/**
 * @brief Tell the #qsched that a task has completed.
 *
 * @param s Pointer to the #qsched.
 * @param t Pointer to the completed #task.
 */

void qsched_done(struct qsched *s, struct task *t) {

  int k;
  struct task *t2;

  TIMER_TIC

  /* Set the task stats. */
  t->toc = getticks();
  if (!(s->flags & qsched_flag_norecost))
    t->cost = t->toc - t->tic;

  /* Release this task's locks. */
  for (k = 0; k < t->nr_locks; k++)
    qsched_unlockres(s, t->locks[k]);

  /* Loop over the task's unlocks... */
  for (k = 0; k < t->nr_unlocks; k++) {

    /* Get a grip on the unlocked task. */
    t2 = &s->tasks[t->unlocks[k]];

    /* Is the unlocked task ready to run? */
    if (atomic_dec(&t2->wait) == 1 && !(t2->flags & task_flag_skip))
      qsched_enqueue(s, t2);
  }

  /* Decrease the number of tasks in this space. */
  atomic_dec(&s->waiting);

/* Ring a bell? */
#ifdef HAVE_PTHREAD
  if (s->flags & qsched_flag_yield) {
    pthread_mutex_lock(&s->mutex);
    pthread_cond_broadcast(&s->cond);
    pthread_mutex_unlock(&s->mutex);
  }
#endif

  /* Careful, this may pick up duplicate timers if virtual
     tasks are used. */
  TIMER_TOC(s, qsched_timer_done);
}

/**
 * @brief Lock a resource and hold its parents.
 *
 * @param s Pointer to the #qsched.
 * @param rid The ID of the resource to lock.
 *
 * @return @c 1 if the resource could be locked, @c 0 otherwise.
 */

int qsched_lockres(struct qsched *s, int rid) {

  int finger, finger2;

  /* Try to lock the root-level resource. */
  if (s->res[rid].hold || lock_trylock(&s->res[rid].lock))
    return 0;

  /* Did the resource get held in the meantime? */
  if (s->res[rid].hold) {
    lock_unlock_blind(&s->res[rid].lock);
    return 0;
  }

  /* Follow parents and increase their hold counter, but fail
     if any are locked. */
  for (finger = s->res[rid].parent; finger != qsched_res_none;
       finger = s->res[finger].parent) {
    if (lock_trylock(&s->res[finger].lock))
      break;
    atomic_inc(&s->res[finger].hold);
    lock_unlock_blind(&s->res[finger].lock);
  }

  /* Did we fail on the way up? */
  if (finger != qsched_res_none) {

    /* Unlock the resource. */
    lock_unlock_blind(&s->res[rid].lock);

    /* Go back up the tree and undo the holds. */
    for (finger2 = s->res[rid].parent; finger2 != finger;
         finger2 = s->res[finger2].parent)
      atomic_dec(&s->res[finger2].hold);

    /* Fail. */
    return 0;

  }

  /* Otherwise, all went well. */
  else
    return 1;
}

/**
 * @brief Unlock a resource and un-hold its parents.
 *
 * @param s Pointer to the #qsched.
 * @param rid The ID of the resource to lock.
 */

void qsched_unlockres(struct qsched *s, int rid) {

  int finger;

  /* Unlock the resource. */
  lock_unlock_blind(&s->res[rid].lock);

  /* Go back up the tree and undo the holds. */
  for (finger = s->res[rid].parent; finger != qsched_res_none;
       finger = s->res[finger].parent)
    atomic_dec(&s->res[finger].hold);
}

/**
 * @brief Try to get all the locks for a task.
 *
 * @param s Pointer to the #qsched.
 * @param tid The ID of the #task to lock.
 *
 * @return @c 1 if the resources could be locked, @c 0 otherwise.
 */

int qsched_locktask(struct qsched *s, int tid) {

  int k;
  struct task *t;

  TIMER_TIC

  /* Get a pointer on the task. */
  t = &s->tasks[tid];

  /* Try to lock all the task's locks. */
  for (k = 0; k < t->nr_locks; k++)
    if (qsched_lockres(s, t->locks[k]) == 0)
      break;

  /* If I didn't get all the locks... */
  if (k < t->nr_locks) {

    /* Unroll the locks I got. */
    for (k -= 1; k >= 0; k--)
      qsched_unlockres(s, t->locks[k]);

    /* Fail. */
    TIMER_TOC(s, qsched_timer_lock);
    return 0;

  }

  /* Otherwise, all went well. */
  else {
    TIMER_TOC(s, qsched_timer_lock);
    return 1;
  }
}

/**
 * @brief Unlock the resources associated with a task.
 *
 * @param s Pointer to the #qsched.
 * @param tid The ID of the #task to unlock.
 */

void qsched_unlocktask(struct qsched *s, int tid) {

  int k;
  struct task *t;

  TIMER_TIC

  /* Get a pointer on the task. */
  t = &s->tasks[tid];

  /* Unlock the used resources. */
  for (k = 0; k < t->nr_locks; k++)
    qsched_unlockres(s, t->locks[k]);

  TIMER_TOC(s, qsched_timer_lock);
}

/**
 * @brief Get a task from the #qsched.
 *
 * @param s Pointer to the #qsched.
 * @param qid The queue to use.
 *
 * @return A pointer to a task object.
 *
 * Note that the #qsched has to have been prepared with #qsched_prepare
 * before any tasks can be extracted. Adding dependencies or locks
 * will require the #qsched to be re-prepared.
 */

struct task *qsched_gettask(struct qsched *s, int qid) {

  int naq, k, tid, qids[s->nr_queues];
  struct task *t;
  unsigned int seed = qid;

  TIMER_TIC

  /* Check if the sched is ok. */
  if (s->flags & qsched_flag_dirty || !(s->flags & qsched_flag_ready))
    quicksched_error("Calling gettask with dirty or unprepared sched.");

  /* Check if the queue ID is ok. */
  if (qid < 0 || qid >= s->nr_queues)
    quicksched_error("Invalid queue ID.");

  /* Main loop. */
  while (s->waiting) {

    /* Try to get a task from my own queue. */
    {
      TIMER_TIC
      tid = queue_get(&s->queues[qid], s, 1);
      TIMER_TOC(s, qsched_timer_queue);
      if (tid < 0) {

        /* Otherwise, hit the other queues. */
        for (naq = 0, k = 0; k < s->nr_queues; k++)
          if (k != qid && s->queues[k].count > 0)
            qids[naq++] = k;
        while (naq > 0) {
          k = rand_r(&seed) % naq;
          TIMER_TIC2
          tid = queue_get(&s->queues[qids[k]], s, 0);
          TIMER_TOC(s, qsched_timer_queue);
          if (tid < 0)
            qids[k] = qids[--naq];
          else
            break;
        }
      }
    }

    /* Bail if a valid task ID was returned. */
    if (tid >= 0) {

      /* Get a pointer to the task. */
      t = &s->tasks[tid];

      /* Own the resources. */
      if (!(s->flags & qsched_flag_noreown)) {
        for (k = 0; k < t->nr_locks; k++)
          s->res[t->locks[k]].owner = qid;
        for (k = 0; k < t->nr_uses; k++)
          s->res[t->uses[k]].owner = qid;
      }

      /* Set some stats data. */
      t->tic = getticks();
      t->qid = qid;

      /* Return the task. */
      TIMER_TOC(s, qsched_timer_gettask);
      return t;

    }

/* Otherwise, take a nap? */
#ifdef HAVE_PTHREAD
    else if (s->flags & qsched_flag_yield) {
      TIMER_TOC(s, qsched_timer_gettask);
      pthread_mutex_lock(&s->mutex);
      if (s->waiting)
        pthread_cond_wait(&s->cond, &s->mutex);
      pthread_mutex_unlock(&s->mutex);
      TIMER_TIC2
    }
#endif
  }

  /* Return empty-handed. No toc here as we don't want to
     count the final wait when all tasks have been executed. */
  return NULL;
}

/**
 * @brief Sort the data according to the given indices.
 *
 * @param data The data to be sorted
 * @param ind The indices with respect to which the data are sorted.
 * @param N The number of entries
 * @param min Lowest index.
 * @param max highest index.
 */

void qsched_sort(int *restrict data, int *restrict ind, int N, int min,
                 int max) {
  int *new_data;
  int *new_ind;
  int i;
  if (N <= 0)
    return;
  new_data = (int *)malloc(sizeof(int) * N);
  if (new_data == NULL)
    quicksched_error("Failed to allocate new_data");
  new_ind = (int *)malloc(sizeof(int) * N);
  if (new_ind == NULL)
    quicksched_error("Failed to allocate new_ind");

  /*Create buckets of size ? - Ideally <16 elements per bucket. Use max-min / N
   * * 10 ? Should give average of 10 elements per bucket */
  int bucketsize = 1;

  /* To find bucket do ind-min / b and it goes in that bucket.*/
  int num_buckets = (max - min) / bucketsize + 1;
  int *bucket_inds = (int *)malloc(sizeof(int) * num_buckets);
  if (bucket_inds == NULL)
    quicksched_error("Failed to allocate bucket_inds");
  memset(bucket_inds, 0, sizeof(int) * num_buckets);
  for (i = 0; i < N; i++) {
    bucket_inds[(ind[i] - min)]++;
  }
  for (i = 1; i < num_buckets; i++) {
    bucket_inds[i] = bucket_inds[i] + bucket_inds[i - 1];
  }
  /* bucket_inds[i] contains the starting position for the i+1th bucket*/
  for (i = num_buckets - 1; i > 0; i--) {
    bucket_inds[i] = bucket_inds[i - 1];
  }
  bucket_inds[0] = 0;

  for (i = 0; i < N; i++) {
    int z = (ind[i] - min);
    new_data[bucket_inds[z]] = data[i];
    new_ind[bucket_inds[z]++] = ind[i];
  }

  /* Copy data back to data and ind and deallocate everything!*/
  memcpy(data, new_data, N * sizeof(int));
  memcpy(ind, new_ind, N * sizeof(int));
  free(new_data);
  free(new_ind);
  free(bucket_inds);
}

/**
 * @brief Sort the data according to the given indices.
 *
 * @param data The data to be sorted
 * @param ind The indices with respect to which the data are sorted.
 * @param N The number of entries
 * @param min Lowest index.
 * @param max highest index.
 *
 * This function calls itself recursively.
 */

void qsched_quicksort(int *restrict data, int *restrict ind, int N, int min,
                      int max) {

  int pivot = (min + max) / 2;
  int i = 0, j = N - 1;
  int temp_i, temp_d;

  /* If N is small enough, just do insert sort. */
  if (N < 16) {

    for (i = 1; i < N; i++)
      if (ind[i] < ind[i - 1]) {
        temp_i = ind[i];
        temp_d = data[i];
        for (j = i; j > 0 && ind[j - 1] > temp_i; j--) {
          ind[j] = ind[j - 1];
          data[j] = data[j - 1];
        }
        ind[j] = temp_i;
        data[j] = temp_d;
      }

  }

  /* Otherwise, recurse with Quicksort. */
  else {

    /* One pass of quicksort. */
    while (i < j) {
      while (i < N && ind[i] <= pivot)
        i++;
      while (j >= 0 && ind[j] > pivot)
        j--;
      if (i < j) {
        temp_i = ind[i];
        ind[i] = ind[j];
        ind[j] = temp_i;
        temp_d = data[i];
        data[i] = data[j];
        data[j] = temp_d;
      }
    }

    /* Recurse in parallel? */
    if (N > 100) {

      /* Recurse on the left? */
      if (j > 0 && pivot > min) {
#pragma omp task untied
        qsched_quicksort(data, ind, j + 1, min, pivot);
      }

      /* Recurse on the right? */
      if (i < N && pivot + 1 < max) {
#pragma omp task untied
        qsched_quicksort(&data[i], &ind[i], N - i, pivot + 1, max);
      }

    } else {

      /* Recurse on the left? */
      if (j > 0 && pivot > min)
        qsched_quicksort(data, ind, j + 1, min, pivot);

      /* Recurse on the right? */
      if (i < N && pivot + 1 < max)
        qsched_quicksort(&data[i], &ind[i], N - i, pivot + 1, max);
    }
  }
}

/**
 * @brief Prepare a #qsched for execution.
 *
 * @param s Pointer to the #qsched.
 */

void qsched_prepare(struct qsched *s) {

  int j, k, count;
  struct task *t, *tasks;

  TIMER_TIC

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Get a pointer to the tasks, set the count. */
  tasks = s->tasks;
  count = s->count;
  /* If the sched is dirty... */
  if (s->flags & qsched_flag_dirty) {

/* Do the sorts in parallel, if possible. */
#pragma omp parallel
    {

/* Sort the unlocks. */
#pragma omp single nowait
      qsched_sort(s->deps, s->deps_key, s->count_deps, 0, count - 1);

/* Sort the locks. */
#pragma omp single nowait
      qsched_sort(s->locks, s->locks_key, s->count_locks, 0, count - 1);

/* Sort the uses. */
#pragma omp single nowait
      qsched_sort(s->uses, s->uses_key, s->count_uses, 0, count - 1);
    }
    /* Run throught the tasks and link the locks and unlocks. */
    tasks[0].unlocks = s->deps;
    tasks[0].locks = s->locks;
    tasks[0].uses = s->uses;
    for (k = 1; k < count; k++) {
      tasks[k].unlocks = &tasks[k - 1].unlocks[tasks[k - 1].nr_unlocks];
      tasks[k].locks = &tasks[k - 1].locks[tasks[k - 1].nr_locks];
      tasks[k].uses = &tasks[k - 1].uses[tasks[k - 1].nr_uses];
    }

    /* All cleaned-up now! */
    s->flags &= ~qsched_flag_dirty;
  }

  /* Init the queues. */
  for (k = 0; k < s->nr_queues; k++)
    queue_init(&s->queues[k], count);

  /* Run through the tasks and set the waits... */
  for (k = 0; k < count; k++) {
    t = &tasks[k];
    if (!(t->flags & task_flag_skip))
      for (j = 0; j < t->nr_unlocks; j++)
        tasks[t->unlocks[j]].wait += 1;
  }

  /* Sort the tasks topologically. */
  int *tid = (int *)malloc(sizeof(int) * count);
  for (j = 0, k = 0; k < count; k++)
    if (tasks[k].wait == 0) {
      tid[j] = k;
      j += 1;
    }
  int ready = j;
  for (k = 0; k < j; k++) {
    t = &tasks[tid[k]];
    for (int kk = 0; kk < t->nr_unlocks; kk++)
      if ((tasks[t->unlocks[kk]].wait -= 1) == 0) {
        tid[j] = t->unlocks[kk];
        j += 1;
      }
  }
  if (k < count)
    quicksched_error("Circular dependencies detected.");

  /* Run through the topologically sorted tasks backwards and
     set their weights, re-setting the waits while we're at it. */
  for (k = count - 1; k >= 0; k--) {
    long long int maxweight = 0;
    t = &tasks[tid[k]];
    for (j = 0; j < t->nr_unlocks; j++) {
      tasks[t->unlocks[j]].wait += 1;
      if (tasks[t->unlocks[j]].weight > maxweight)
        maxweight = tasks[t->unlocks[j]].weight;
    }
    t->weight = t->cost + maxweight;
  }

  /* Run through the tasks and enqueue the non-waiting ones. */
  for (k = 0; k < ready; k++) {
    t = &tasks[tid[k]];
    if (t->wait == 0 && !(t->flags & task_flag_skip))
      qsched_enqueue(s, t);
  }

  /* Clean up. */
  free(tid);

  /* Set the number of waiting tasks. */
  s->waiting = count;

  /* Set the ready flag. */
  s->flags |= qsched_flag_ready;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);

  TIMER_TOC(s, qsched_timer_prepare);
}

/**
 * @brief Add a new resource to the #qsched.
 *
 * @param s Pointer to the #qsched.
 * @param parent ID of the parent resource or #qsched_res_none if none.
 * @param owner ID of the ower
 *
 * @return The ID of the new shared resource.
 */

int qsched_addres(struct qsched *s, int owner, int parent) {

  struct res *res_new;
  int id;

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Do the deps need to be re-allocated? */
  if (s->count_res == s->size_res) {

    /* Scale the res list size. */
    s->size_res *= qsched_stretch;

    /* Allocate a new task list. */
    if ((res_new = malloc(sizeof(struct res) * s->size_res)) == NULL)
      quicksched_error("Failed to allocate new res lists.");

    /* Copy the res and owners over to the new list. */
    memcpy(res_new, s->res, sizeof(struct res) * s->count_res);

    /* Free the old res lists. */
    free(s->res);

    /* Set the new res lists. */
    s->res = res_new;
  }

  /* Increase the res counter. */
  id = s->count_res;
  s->count_res += 1;

  /* Init the resource. */
  lock_init(&s->res[id].lock);
  s->res[id].hold = 0;
  s->res[id].owner = owner;
  s->res[id].parent = parent;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);

  /* Return the res ID. */
  return id;
}

/**
 * @brief Add a resource requirement to a task.
 *
 * @param s Pointer to the #qsched.
 * @param t ID of the task.
 * @param res ID of the resource.
 */

void qsched_addlock(struct qsched *s, int t, int res) {

  void *temp1, *temp2;

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Do the deps need to be re-allocated? */
  if (s->count_locks == s->size_locks) {

    /* Scale the locks list size. */
    s->size_locks *= qsched_stretch;

    /* Allocate a new task list. */
    if ((temp1 = malloc(sizeof(int) * s->size_locks)) == NULL ||
        (temp2 = malloc(sizeof(int) * s->size_locks)) == NULL)
      quicksched_error("Failed to allocate new locks lists.");

    /* Copy the locks and keys over to the new list. */
    memcpy(temp1, s->locks, sizeof(int) * s->count_locks);
    memcpy(temp2, s->locks_key, sizeof(int) * s->count_locks);

    /* Free the old locks lists. */
    free(s->locks);
    free(s->locks_key);

    /* Set the new locks lists. */
    s->locks = (int *)temp1;
    s->locks_key = (int *)temp2;
  }

  /* Add the new dependency. */
  s->locks[s->count_locks] = res;
  s->locks_key[s->count_locks] = t;
  s->tasks[t].nr_locks += 1;

  /* Increase the locks counter. */
  s->count_locks += 1;

  /* The sched is now dirty. */
  s->flags |= qsched_flag_dirty;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Add a resource use to a task.
 *
 * @param s Pointer to the #qsched.
 * @param t ID of the task.
 * @param res ID of the resource.
 */

void qsched_adduse(struct qsched *s, int t, int res) {

  void *temp1, *temp2;

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Do the deps need to be re-allocated? */
  if (s->count_uses == s->size_uses) {

    /* Scale the uses list size. */
    s->size_uses *= qsched_stretch;

    /* Allocate a new task list. */
    if ((temp1 = malloc(sizeof(int) * s->size_uses)) == NULL ||
        (temp2 = malloc(sizeof(int) * s->size_uses)) == NULL)
      quicksched_error("Failed to allocate new uses lists.");

    /* Copy the uses and keys over to the new list. */
    memcpy(temp1, s->uses, sizeof(int) * s->count_uses);
    memcpy(temp2, s->uses_key, sizeof(int) * s->count_uses);

    /* Free the old uses lists. */
    free(s->uses);
    free(s->uses_key);

    /* Set the new uses lists. */
    s->uses = (int *)temp1;
    s->uses_key = (int *)temp2;
  }

  /* Add the new dependency. */
  s->uses[s->count_uses] = res;
  s->uses_key[s->count_uses] = t;
  s->tasks[t].nr_uses += 1;

  /* Increase the uses counter. */
  s->count_uses += 1;

  /* The sched is now dirty. */
  s->flags |= qsched_flag_dirty;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Add a task dependency.
 *
 * @param s Pointer to the #qsched.
 * @param ta ID of the unlocking task.
 * @param tb ID of the unlocked task.
 *
 * A dependency is added such that @c tb depends on @c ta.
 */

void qsched_addunlock(struct qsched *s, int ta, int tb) {

  void *temp1, *temp2;

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Do the deps need to be re-allocated? */
  if (s->count_deps == s->size_deps) {

    /* Scale the deps list size. */
    s->size_deps *= qsched_stretch;

    /* Allocate a new task list. */
    if ((temp1 = malloc(sizeof(int) * s->size_deps)) == NULL ||
        (temp2 = malloc(sizeof(int) * s->size_deps)) == NULL)
      quicksched_error("Failed to allocate new deps lists.");

    /* Copy the deps and keys over to the new list. */
    memcpy(temp1, s->deps, sizeof(int) * s->count_deps);
    memcpy(temp2, s->deps_key, sizeof(int) * s->count_deps);

    /* Free the old deps lists. */
    free(s->deps);
    free(s->deps_key);

    /* Set the new deps lists. */
    s->deps = (int *)temp1;
    s->deps_key = (int *)temp2;
  }

  /* Add the new dependency. */
  s->deps[s->count_deps] = tb;
  s->deps_key[s->count_deps] = ta;
  s->tasks[ta].nr_unlocks += 1;

  /* Increase the deps counter. */
  s->count_deps += 1;

  /* The sched is now dirty. */
  s->flags |= qsched_flag_dirty;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Add a new task to the #qsched.
 *
 * @param s Pointer to the #qsched
 * @param type Task type.
 * @param flags Task flags.
 * @param data Pointer to the task data.
 * @param data_size Size, in bytes, of the task data.
 * @param cost Approximate cost for this task.
 */

int qsched_addtask(struct qsched *s, int type, unsigned int flags, void *data,
                   int data_size, int cost) {

  void *temp;
  struct task *t;
  int id, data_size2;

  /* Lock the sched. */
  lock_lock(&s->lock);

  /* Do the tasks need to be re-allocated? */
  if (s->count == s->size) {

    /* Scale the task list size. */
    s->size *= qsched_stretch;

    /* Allocate a new task list. */
    if ((temp = malloc(sizeof(struct task) * s->size)) == NULL)
      quicksched_error("Failed to allocate new task list.");

    /* Copy the tasks over to the new list. */
    memcpy(temp, s->tasks, sizeof(struct task) * s->count);

    /* Free the old task list. */
    free(s->tasks);

    /* Set the new task list. */
    s->tasks = (struct task *)temp;
  }

  /* Round-up the data size. */
  data_size2 = (data_size + (qsched_data_round - 1)) & ~(qsched_data_round - 1);

  /* Do the task data need to be re-allocated? */
  if (s->count_data + data_size2 > s->size_data) {

    /* Scale the task list size. */
    s->size_data *= qsched_stretch;

    /* Allocate a new task list. */
    if ((temp = malloc(s->size_data)) == NULL)
      quicksched_error("Failed to allocate new task list.");

    /* Copy the tasks over to the new list. */
    memcpy(temp, s->data, s->count_data);

    /* Free the old task list. */
    free(s->data);

    /* Set the new task list. */
    s->data = temp;
  }

  /* Store the new task ID. */
  id = s->count;

  /* Init the new task. */
  t = &s->tasks[id];
  t->type = type;
  t->flags = flags;
  t->cost = cost;
  t->wait = 0;
  t->nr_conflicts = 0;
  t->nr_unlocks = 0;
  t->nr_locks = 0;
  t->nr_uses = 0;

  /* Add a relative pointer to the data. */
  memcpy(&s->data[s->count_data], data, data_size);
  t->data = &s->data[s->count_data] - s->data;
  s->count_data += data_size2;

  /* Increase the task counter. */
  s->count += 1;

  /* Unlock the sched. */
  lock_unlock_blind(&s->lock);

  /* Return the task ID. */
  return id;
}

/**
 * @brief Clean up a #qsched, free all associated memory.
 *
 * @param s Pointer to the #qsched.
 */

void qsched_free(struct qsched *s) {

  int k;

  /* Clear all the buffers if allocated. */
  if (s->tasks != NULL) {
    free(s->tasks);
    s->tasks = NULL;
  }
  if (s->deps != NULL) {
    free(s->deps);
    s->deps = NULL;
  }
  if (s->deps_key != NULL) {
    free(s->deps_key);
    s->deps_key = NULL;
  }
  if (s->locks != NULL) {
    free(s->locks);
    s->locks = NULL;
  }
  if (s->locks_key != NULL) {
    free(s->locks_key);
    s->locks_key = NULL;
  }
  if (s->uses != NULL) {
    free(s->uses);
    s->uses = NULL;
  }
  if (s->uses_key != NULL) {
    free(s->uses_key);
    s->uses_key = NULL;
  }
  if (s->res != NULL) {
    free((void *)s->res);
    s->res = NULL;
  }
  if (s->data != NULL) {
    free(s->data);
    s->data = NULL;
  }

  /* Loop over the queues and free them too. */
  for (k = 0; k < s->nr_queues; k++)
    queue_free(&s->queues[k]);
  free(s->queues);
  s->queues = NULL;

/* Destroy the mutex and condition. */
#ifdef HAVE_PTHREAD
  if (s->flags & qsched_flag_pthread) {

    /* Start all the threads on an empty function, to kill them. */
    s->fun = NULL;
    s->barrier_count = s->runners_count;
    s->barrier_launchcount = s->runners_count;
    if (pthread_mutex_unlock(&s->barrier_mutex) != 0 ||
        pthread_cond_broadcast(&s->barrier_cond) != 0)
      quicksched_error("Failed to open the barrier.");

    /* Wait for each thread to have terminated. */
    for (k = 0; k < s->runners_count; k++)
      if (pthread_join(s->runners[k].thread, NULL) != 0)
        quicksched_error("Failed to join on thread %i.", k);

    /* Clean up the mutexes and barriers. */
    if (pthread_cond_destroy(&s->cond) != 0 ||
        pthread_mutex_destroy(&s->mutex) != 0)
      quicksched_error("Error destroying pthread cond/mutex pair.");
    if (pthread_mutex_destroy(&s->barrier_mutex) != 0 ||
        pthread_cond_destroy(&s->barrier_cond) != 0)
      quicksched_error("Error destroying pthread barrier cond/mutex pair.");
    free(s->runners);
    s->runners_size = 0;
    s->runners_count = 0;
  }
#endif

  /* Clear the flags. */
  s->flags = qsched_flag_none;
}

/**
 * @brief Initialize the given #qsched object.
 *
 * @param s Pointer to a #qsched object.
 * @param nr_queues The number of queues in the #qsched.
 * @param flags Flags specifying the behaviour of this #qsched.
 *
 * Initializes the given #qsched with the given number of queues.
 */

void qsched_init(struct qsched *s, int nr_queues, int flags) {

  /* Set the flags to begin with. */
  s->flags = flags;

  /* Allocate and clear the queues (will init when sched is
     finalized. */
  if ((s->queues = (struct queue *)malloc(sizeof(struct queue) * nr_queues)) ==
      NULL)
    quicksched_error("Failed to allocate memory for queues.");
  bzero(s->queues, sizeof(struct queue) * nr_queues);
  s->nr_queues = nr_queues;

  /* Allocate the task list. */
  s->size = qsched_size_init;
  if ((s->tasks = (struct task *)malloc(sizeof(struct task) * s->size)) == NULL)
    quicksched_error("Failed to allocate memory for tasks.");
  s->count = 0;

  /* Allocate the initial deps. */
  s->size_deps = qsched_init_depspertask * s->size;
  if ((s->deps = (int *)malloc(sizeof(int) * s->size_deps)) == NULL ||
      (s->deps_key = (int *)malloc(sizeof(int) * s->size_deps)) == NULL)
    quicksched_error("Failed to allocate memory for deps.");
  s->count_deps = 0;

  /* Allocate the initial locks. */
  s->size_locks = qsched_init_lockspertask * s->size;
  if ((s->locks = (int *)malloc(sizeof(int) * s->size_locks)) == NULL ||
      (s->locks_key = (int *)malloc(sizeof(int) * s->size_locks)) == NULL)
    quicksched_error("Failed to allocate memory for locks.");
  s->count_locks = 0;

  /* Allocate the initial res. */
  s->size_res = qsched_init_respertask * s->size;
  if ((s->res = (struct res *)malloc(sizeof(struct res) * s->size_res)) == NULL)
    quicksched_error("Failed to allocate memory for res.");
  s->count_res = 0;

  /* Allocate the initial uses. */
  s->size_uses = qsched_init_usespertask * s->size;
  if ((s->uses = (int *)malloc(sizeof(int) * s->size_uses)) == NULL ||
      (s->uses_key = (int *)malloc(sizeof(int) * s->size_uses)) == NULL)
    quicksched_error("Failed to allocate memory for uses.");
  s->count_uses = 0;

  /* Allocate the initial data. */
  s->size_data = qsched_init_datapertask * s->size;
  if ((s->data = malloc(s->size_data)) == NULL)
    quicksched_error("Failed to allocate memory for data.");
  s->count_data = 0;

/* Init the pthread stuff. */
#ifdef HAVE_PTHREAD
  if (flags & qsched_flag_pthread) {
    if (pthread_cond_init(&s->cond, NULL) != 0 ||
        pthread_mutex_init(&s->mutex, NULL) != 0)
      quicksched_error("Error initializing yield cond/mutex pair.");
    if (pthread_cond_init(&s->barrier_cond, NULL) != 0 ||
        pthread_mutex_init(&s->barrier_mutex, NULL) != 0)
      quicksched_error("Error initializing barrier cond/mutex pair.");
    s->runners_count = 0;
    s->runners_size = qsched_init_runners;
    if ((s->runners = malloc(sizeof(struct qsched_pthread_runner) *
                             s->runners_size)) == NULL)
      quicksched_error("Failed to allocate runners.");
    s->barrier_running = 0;
    s->barrier_count = 0;
    s->barrier_launchcount = 0;
    if (pthread_mutex_lock(&s->barrier_mutex) != 0)
      quicksched_error("Failed to lock barrier mutex.");
  }
#endif

/* Clear the timers. */
#ifdef TIMERS
  bzero(s->timers, sizeof(ticks) * qsched_timer_count);
#endif

  /* Init the sched lock. */
  lock_init(&s->lock);
}
