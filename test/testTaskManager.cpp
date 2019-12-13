/**
 * @file testTaskManager.cpp
 *
 * @brief Unit test for the TaskManager class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "TaskManager.hpp"

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
 * @brief Unit test for the TaskManager class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, 0);

  task_manager.add_composition(1);
  task_manager.add_size(1.e-6);
  task_manager.add_wavelength(1.e-4);

  QuickSched quicksched(4, true, "test_TaskManager.log");

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  std::vector<TMatrixResource *> tmatrices;
  std::vector<InteractionVariables *> interaction_variables;
  task_manager.generate_tasks(quicksched, tasks, resources, results,
                              space_manager, tmatrices, interaction_variables);

  quicksched.execute_tasks();

  std::ofstream taskfile("test_taskmanager_tasks.txt");
  taskfile << "# thread\tstart\tend\ttype\ttask id\n";
  print_vector(tasks, quicksched, taskfile);
  std::ofstream typefile("test_taskmanager_types.txt");
  typefile << "# type\tlabel\n";
  quicksched.print_type_dict(typefile);

  clear_vector(tasks);
  clear_vector(resources);
  clear_vector(results);
  delete space_manager;
  clear_vector(tmatrices);
  clear_vector(interaction_variables);

  return 0;
}
