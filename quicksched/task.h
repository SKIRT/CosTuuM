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

/* Task flags. */
#define task_flag_none 0
#define task_flag_skip 1
#define task_flag_virtual 2

/* The task data structure. */
struct task {

  /* Task type. */
  int type;

  /* Task flags. */
  unsigned int flags;

  /* Task payload offset. */
  int data;

  /* Task wait counter. */
  int wait;

  /* Number of potential conflicts. */
  int nr_conflicts;

  /* Task unlocks. */
  int *unlocks, nr_unlocks;

  /* Task locks. */
  int *locks, nr_locks;

  /* Task uses. */
  int *uses, nr_uses;

  /* Timers and other info. */
  ticks tic, toc;
  int qid;

  /* Task weight for queue selection. */
  int cost, weight;
};
