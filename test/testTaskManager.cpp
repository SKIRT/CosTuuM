/**
 * @file testTaskManager.cpp
 *
 * @brief Unit test for the TaskManager class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "TaskManager.hpp"

#include <fstream>
#include <string.h>

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
    // make sure we used all elements in the vectors
    ctm_assert(vec[i] != nullptr);
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
 * @brief Write the given integer to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 * @tparam T Integer type.
 */
template <typename T>
inline void int_to_file(std::ofstream &ofile, const T value) {
  const uint_fast64_t intval = value;
  ofile.write(reinterpret_cast<const char *>(&intval), sizeof(intval));
}

/**
 * @brief Write the given string to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 */
inline void string_to_file(std::ofstream &ofile, const std::string value) {
  char strval[] = "        ";
  memcpy(strval, value.c_str(), value.size());
  ofile.write(strval, 8);
}

/**
 * @brief Write the given array to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 */
inline void array_to_file(std::ofstream &ofile,
                          const std::vector<float_type> &value) {
  for (uint_fast32_t i = 0; i < value.size(); ++i) {
    const double ival = value[i];
    ofile.write(reinterpret_cast<const char *>(&ival), sizeof(ival));
  }
}

/**
 * @brief Write the given float to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 * @tparam T Floating point type.
 */
template <typename T>
inline void float_to_file(std::ofstream &ofile, const T &value) {
  const double fltval = value;
  ofile.write(reinterpret_cast<const char *>(&fltval), sizeof(fltval));
}

/**
 * @brief Unit test for the TaskManager class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, 1, 0);

  std::vector<float_type> sizes(2);
  sizes[0] = 1.e-7;
  sizes[1] = 1.e-5;

  std::vector<float_type> wavelengths(3);
  wavelengths[0] = 1.e-4;
  wavelengths[1] = 1.e-3;
  wavelengths[2] = 1.e-2;

  task_manager.add_composition(DUSTGRAINTYPE_SILICON);
  for (uint_fast32_t i = 0; i < sizes.size(); ++i) {
    task_manager.add_size(sizes[i]);
  }
  for (uint_fast32_t i = 0; i < wavelengths.size(); ++i) {
    task_manager.add_wavelength(wavelengths[i]);
  }

  QuickSched quicksched(4, true, "test_TaskManager.log");

  const uint_fast32_t ntheta = 10;
  std::vector<float_type> thetas(ntheta);
  for (uint_fast32_t i = 0; i < ntheta; ++i) {
    thetas[i] = (i + 0.5) * M_PI / ntheta;
  }

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_tasks(thetas, 20, quicksched, tasks, resources, results,
                              space_manager);

  quicksched.execute_tasks();

  std::ofstream taskfile("test_taskmanager_tasks.txt");
  taskfile << "# thread\tstart\tend\ttype\ttask id\n";
  print_vector(tasks, quicksched, taskfile);
  std::ofstream typefile("test_taskmanager_types.txt");
  typefile << "# type\tlabel\n";
  quicksched.print_type_dict(typefile);

  // output some example results along each axis
  {
    std::ofstream ofile("test_taskmanager_size.txt");
    ofile << "# size\tQabs\tQabspol\n";
    for (uint_fast32_t isize = 0; isize < sizes.size(); ++isize) {
      const uint_fast32_t result_index = isize * wavelengths.size();
      const AbsorptionCoefficientResult &result =
          *static_cast<AbsorptionCoefficientResult *>(results[result_index]);
      ofile << sizes[isize] << "\t" << result.get_Qabs(0) << "\t"
            << result.get_Qabspol(0) << "\n";
    }
  }
  {
    std::ofstream ofile("test_taskmanager_wavelength.txt");
    ofile << "# wavelength\tQabs\tQabspol\n";
    for (uint_fast32_t ilambda = 0; ilambda < wavelengths.size(); ++ilambda) {
      const uint_fast32_t result_index = ilambda;
      const AbsorptionCoefficientResult &result =
          *static_cast<AbsorptionCoefficientResult *>(results[result_index]);
      ofile << wavelengths[ilambda] << "\t" << result.get_Qabs(0) << "\t"
            << result.get_Qabspol(0) << "\n";
    }
  }
  {
    std::ofstream ofile("test_taskmanager_theta.txt");
    ofile << "# theta\tQabs\tQabspol\n";
    const AbsorptionCoefficientResult &result =
        *static_cast<AbsorptionCoefficientResult *>(results[0]);
    for (uint_fast32_t itheta = 0; itheta < thetas.size(); ++itheta) {
      ofile << thetas[itheta] << "\t" << result.get_Qabs(itheta) << "\t"
            << result.get_Qabspol(itheta) << "\n";
    }
  }

  {
    std::ofstream ofile("test_table.stab", std::ios::binary);
    const std::string skirt_tag("SKIRT X\n");
    ofile.write(skirt_tag.c_str(), skirt_tag.size());
    int_to_file(ofile, 0x010203040A0BFEFF);

    int_to_file(ofile, 3);
    string_to_file(ofile, "a");
    string_to_file(ofile, "lambda");
    string_to_file(ofile, "theta");
    string_to_file(ofile, "m");
    string_to_file(ofile, "m");
    string_to_file(ofile, "rad");
    string_to_file(ofile, "log");
    string_to_file(ofile, "log");
    string_to_file(ofile, "lin");

    int_to_file(ofile, sizes.size());
    array_to_file(ofile, sizes);
    int_to_file(ofile, wavelengths.size());
    array_to_file(ofile, wavelengths);
    int_to_file(ofile, thetas.size());
    array_to_file(ofile, thetas);

    int_to_file(ofile, 2);
    string_to_file(ofile, "Qabs");
    string_to_file(ofile, "Qabspol");
    string_to_file(ofile, "1");
    string_to_file(ofile, "1");
    string_to_file(ofile, "log");
    string_to_file(ofile, "lin");
    for (uint_fast32_t itheta = 0; itheta < thetas.size(); ++itheta) {
      for (uint_fast32_t ilambda = 0; ilambda < wavelengths.size(); ++ilambda) {
        for (uint_fast32_t isize = 0; isize < sizes.size(); ++isize) {
          const uint_fast32_t result_index =
              isize * wavelengths.size() + ilambda;
          const AbsorptionCoefficientResult &result =
              *static_cast<AbsorptionCoefficientResult *>(
                  results[result_index]);
          assert_condition(sizes[isize] == result.get_size());
          assert_condition(wavelengths[ilambda] == result.get_wavelength());
          float_to_file(ofile, result.get_Qabs(itheta));
          float_to_file(ofile, result.get_Qabspol(itheta));
        }
      }
    }

    const std::string trail_tag("STABEND\n");
    ofile.write(trail_tag.c_str(), trail_tag.size());
  }

  clear_vector(tasks);
  clear_vector(resources);
  clear_vector(results);
  delete space_manager;

  return 0;
}
