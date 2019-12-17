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

  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, 2, 0);

  task_manager.add_composition(DUSTGRAINTYPE_SILICON);
  task_manager.add_size(1.e-7);
  //  task_manager.add_size(1.e-5);
  task_manager.add_wavelength(1.e-4);
  //  task_manager.add_wavelength(1.e-3);

  QuickSched quicksched(4, true, "test_TaskManager.log");

  const uint_fast32_t ntheta = 10;
  std::vector<float_type> thetas(ntheta);
  for (uint_fast32_t i = 0; i < ntheta; ++i) {
    thetas[i] = (i + 0.5) * M_PI / ntheta;
  }
  AbsorptionCoefficientGrid grid(ntheta, &thetas[0]);

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_tasks(grid, quicksched, tasks, resources, results,
                              space_manager);

  quicksched.execute_tasks();

  std::ofstream taskfile("test_taskmanager_tasks.txt");
  taskfile << "# thread\tstart\tend\ttype\ttask id\n";
  print_vector(tasks, quicksched, taskfile);
  std::ofstream typefile("test_taskmanager_types.txt");
  typefile << "# type\tlabel\n";
  quicksched.print_type_dict(typefile);

  const AbsorptionCoefficientResult &result =
      *static_cast<AbsorptionCoefficientResult *>(results[0]);
  ctm_warning("composition: %" PRIiFAST32, result.get_composition());
  ctm_warning("size: %g", double(result.get_size()));
  ctm_warning("wavelength: %g", double(result.get_wavelength()));
  ctm_warning("type: %" PRIiFAST32, result.get_type());
  for (uint_fast32_t i = 0; i < ntheta; ++i) {
    ctm_warning("%g %g %g", double(thetas[i]), double(result.get_Qabs(i)),
                double(result.get_Qabspol(i)));
  }

  clear_vector(tasks);
  clear_vector(resources);
  clear_vector(results);
  delete space_manager;

  return 0;
}
