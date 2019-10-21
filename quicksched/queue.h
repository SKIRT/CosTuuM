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

/* Some constants. */
#define queue_grow 2

/* The queue data structure. */
struct queue {

  /* Task indices. */
  int *inds;

  /* Number of tasks waiting in the queue. */
  int count;

  /* Lock to exclusively access this queue. */
  lock_type lock;

  /* Maximum number of tasks in queue. */
  int size;

} __attribute__((aligned(128)));

/* Function prototypes. */
int queue_get(struct queue *q, struct qsched *s, int insist);
void queue_put(struct queue *q, struct qsched *s, int tid);
void queue_init(struct queue *q, int size);
void queue_free(struct queue *q);
