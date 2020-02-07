/**
 * @file CTMModule.cpp
 *
 * @brief CosTuuM Python module.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "PyAlignmentDistribution.hpp"
#include "PyDustProperties.hpp"
#include "PyHelperFunctions.hpp"
#include "PyOrientationDistribution.hpp"
#include "PyShapeDistribution.hpp"
#include "PyTMatrix.hpp"
#include "TaskManager.hpp"

#include <Python.h>
#include <sys/resource.h>
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "structmember.h"

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/// CosTuuM module

/**
 * @brief Get the equal volume radius for the particle with the given equal area
 * radius and axis ratio.
 *
 * Required arguments are:
 *  - radius: Equal area radius (in input length units).
 *  - axis_ratio: Axis ratio.
 *
 * @param self Module object being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Equal volume radius, in same length units as input.
 */
static PyObject *get_equal_volume_radius(PyObject *self, PyObject *args,
                                         PyObject *kwargs) {

  // required arguments
  float_type radius;
  float_type axis_ratio;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("radius"), strdup("axis_ratio"), nullptr};

  // allocate temporary variables to store double precision arguments
  double radius_d, axis_ratio_d;
  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd", kwlist, &radius_d,
                                   &axis_ratio_d)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }
  // unpack double precision arguments
  radius = radius_d;
  axis_ratio = axis_ratio_d;

  // return the array
  return PyFloat_FromDouble(static_cast<double>(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          axis_ratio) *
      radius));
}

/**
 * @brief Get a table of absorption coefficients for the given input dust grain
 * sizes, dust grain types, wavelengths and angles, using the given shape and
 * alignment distributions.
 *
 * Required arguments are:
 *  - types: Array of dust grain types.
 *  - sizes: Array of dust grain sizes (in m).
 *  - wavelengths: Array of incoming/outgoing photon wavelengths (in m).
 *  - thetas: Array of zenith angles (in radians).
 *  - shape_distribution: Shape distribution object to use.
 *  - alignment_distribution: Alignment distribution object to use.
 *  - dust_properties: Dust properties object to use.
 *
 * Additional optional arguments are:
 *  - minimum_order: Minimum order of spherical basis function expansion to use.
 *  - maximum_order: Maximum allowed order of spherical basis function
 *  expansion.
 *  - gauss_legendre_factor: Number of Gauss-Legendre quadrature points to use
 *  as a multiplicative factor of the spherical basis function expansion order.
 *  - tolerance: Maximum allowed relative difference between two successive
 *  orders of T-matrix calculations that decides when the calculation is
 *  converged.
 *  - maximum_memory_size: Maximum allowed memory usage of the algorithm (in
 *  bytes).
 *  - number_of_threads: Number of shared-memory threads to use to run the
 *  calculation.
 *  - quicksched_graph_log: Name for the QuickSched graph log file to write (or
 *  None for no graph log).
 *  - quicksched_task_log: Name for the QuickSched task log file to write (or
 *  None for no task log).
 *  - quicksched_task_type_log: Name for the QuickSched task type log file to
 *  write (or None for no task type log).
 *  - memory_log: Name of the memory log file to write (or None for no memory
 *  log).
 *  - do_absorption: Compute absorption coefficients? (default: True)
 *  - do_extinction: Compute extinction coefficients? (default: False)
 *  - do_scattering: Compute scattering matrices? (default: False)
 *  - account_for_scattering: Subtract the directionally averaged scattering
 *  cross sections from the absorption coefficients? (default: False)
 *  - verbose: Display diagnostic output? (default: False)
 *
 * @param self Module object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Nothing.
 */
static PyObject *get_table(PyObject *self, PyObject *args, PyObject *kwargs) {

  // parse arguments //

  // required arguments
  PyArrayObject *input_types;
  PyArrayObject *input_sizes;
  PyArrayObject *input_wavelengths;
  PyArrayObject *input_thetas;
  PyShapeDistribution *shape_distribution_object;
  PyAlignmentDistribution *alignment_distribution_object;
  PyDustProperties *dust_properties_object;

  // optional arguments
  uint_fast32_t input_nmin = 10;
  uint_fast32_t input_nmax = 100;
  uint_fast32_t input_glfac = 2;
  float_type input_tolerance = 1.e-4;
  size_t input_memory_size = 1e10;
  int_fast32_t input_nthread = 4;
  const char *input_graph_log_name = nullptr;
  const char *input_task_log_name = nullptr;
  const char *input_task_type_log_name = nullptr;
  const char *input_memory_log_name = nullptr;
  bool do_absorption = true;
  bool do_extinction = false;
  bool do_scattering = false;
  bool account_for_scattering = false;
  bool verbose = false;

  // list of keywords
  static char *kwlist[] = {strdup("types"),
                           strdup("sizes"),
                           strdup("wavelengths"),
                           strdup("thetas"),
                           strdup("shape_distribution"),
                           strdup("alignment_distribution"),
                           strdup("dust_properties"),
                           strdup("minimum_order"),
                           strdup("maximum_order"),
                           strdup("gauss_legendre_factor"),
                           strdup("tolerance"),
                           strdup("maximum_memory_size"),
                           strdup("number_of_threads"),
                           strdup("quicksched_graph_log"),
                           strdup("quicksched_task_log"),
                           strdup("quicksched_task_type_log"),
                           strdup("memory_log"),
                           strdup("do_absorption"),
                           strdup("do_extinction"),
                           strdup("do_scattering"),
                           strdup("account_for_scattering"),
                           strdup("verbose"),
                           nullptr};

  // placeholders for float_type arguments
  double input_tolerance_d = static_cast<double>(input_tolerance);
  // placeholders for integer arguments
  unsigned int input_nmin_i = input_nmin;
  unsigned int input_nmax_i = input_nmax;
  unsigned int input_glfac_i = input_glfac;
  unsigned long input_memory_size_i = input_memory_size;
  int input_nthread_i = input_nthread;
  // placeholders for bool arguments
  int do_absorption_b = do_absorption;
  int do_extinction_b = do_extinction;
  int do_scattering_b = do_scattering;
  int account_for_scattering_b = account_for_scattering;
  int verbose_b = verbose;
  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "O&O&O&O&OOO|IIIdkizzzzppppp", kwlist,
          PyArray_Converter, &input_types, PyArray_Converter, &input_sizes,
          PyArray_Converter, &input_wavelengths, PyArray_Converter,
          &input_thetas, &shape_distribution_object,
          &alignment_distribution_object, &dust_properties_object,
          &input_nmin_i, &input_nmax_i, &input_glfac_i, &input_tolerance_d,
          &input_memory_size_i, &input_nthread_i, &input_graph_log_name,
          &input_task_log_name, &input_task_type_log_name,
          &input_memory_log_name, &do_absorption_b, &do_extinction_b,
          &do_scattering_b, &account_for_scattering_b, &verbose_b)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }
  // convert float_type arguments
  input_tolerance = input_tolerance_d;
  // convert integer arguments
  input_nmin = input_nmin_i;
  input_nmax = input_nmax_i;
  input_glfac = input_glfac_i;
  input_memory_size = input_memory_size_i;
  input_nthread = input_nthread_i;
  // convert bool arguments
  do_absorption = do_absorption_b;
  do_extinction = do_extinction_b;
  do_scattering = do_scattering_b;
  account_for_scattering = account_for_scattering_b;
  verbose = verbose_b;

  const ShapeDistribution *shape_distribution =
      shape_distribution_object->_shape_distribution;
  const AlignmentDistribution *alignment_distribution =
      alignment_distribution_object->_alignment_distribution;
  const DustProperties *dust_properties =
      dust_properties_object->_dust_properties;
  TaskManager task_manager(input_nmin, input_nmax, input_glfac, input_tolerance,
                           input_memory_size, *shape_distribution,
                           *alignment_distribution, *dust_properties);

  std::vector<int_fast32_t> types =
      unpack_numpy_array<int_fast32_t>(input_types);
  for (uint_fast32_t i = 0; i < types.size(); ++i) {
    task_manager.add_composition(types[i]);
  }
  Py_DECREF(input_types);

  std::vector<float_type> sizes =
      unpack_numpy_array<float_type, double>(input_sizes);
  for (uint_fast32_t i = 0; i < sizes.size(); ++i) {
    task_manager.add_size(sizes[i]);
  }
  Py_DECREF(input_sizes);

  std::vector<float_type> wavelengths =
      unpack_numpy_array<float_type, double>(input_wavelengths);
  for (uint_fast32_t i = 0; i < wavelengths.size(); ++i) {
    task_manager.add_wavelength(wavelengths[i]);
  }
  Py_DECREF(input_wavelengths);

  std::string input_graph_log;
  if (input_graph_log_name != nullptr) {
    input_graph_log = input_graph_log_name;
  }
  QuickSched quicksched(input_nthread, input_graph_log_name != nullptr,
                        input_graph_log, verbose);

  std::vector<float_type> thetas =
      unpack_numpy_array<float_type, double>(input_thetas);
  Py_DECREF(input_thetas);

  std::string input_memory_log;
  if (input_memory_log_name != nullptr) {
    input_memory_log = input_memory_log_name;
  }
  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  ResultKey *result_key = nullptr;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_tasks(
      thetas, 20, quicksched, tasks, resources, result_key, results,
      space_manager, do_extinction, do_absorption, do_scattering,
      account_for_scattering, verbose, input_memory_log_name != nullptr,
      input_memory_log);

  Py_BEGIN_ALLOW_THREADS;
  quicksched.execute_tasks();
  Py_END_ALLOW_THREADS;

  if (verbose) {
    struct rusage resource_usage;
    getrusage(RUSAGE_SELF, &resource_usage);
    size_t memory_usage = static_cast<size_t>(resource_usage.ru_maxrss) *
                          static_cast<size_t>(1024);
    ctm_warning("Actual total memory usage: %s",
                Utilities::human_readable_bytes(memory_usage).c_str());
  }

  if (input_task_log_name) {
    std::ofstream taskfile(input_task_log_name);
    taskfile << "# thread\tstart\tend\ttype\ttask id\n";
    for (uint_fast32_t i = 0; i < tasks.size(); ++i) {
      quicksched.print_task(*tasks[i], taskfile);
    }
  }
  if (input_task_type_log_name) {
    std::ofstream typefile(input_task_type_log_name);
    typefile << "# type\tlabel\n";
    quicksched.print_type_dict(typefile);
  }

  // we are done with the tasks and resources, delete them
  for (uint_fast32_t i = 0; i < tasks.size(); ++i) {
    delete tasks[i];
  }
  for (uint_fast32_t i = 0; i < resources.size(); ++i) {
    delete resources[i];
  }
  delete space_manager;

  PyObject *first_return_value = nullptr;
  if (do_absorption || do_extinction) {
    npy_intp result_size = 0;
    if (do_absorption) {
      result_size += 2;
    }
    if (do_extinction) {
      result_size += 3;
    }
    npy_intp result_dims[5] = {
        static_cast<npy_intp>(result_key->composition_size()),
        static_cast<npy_intp>(result_key->size_size()),
        static_cast<npy_intp>(result_key->wavelength_size()),
        static_cast<npy_intp>(thetas.size()), result_size};
    PyArrayObject *result_array = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(5, result_dims, NPY_DOUBLE));

    for (npy_intp itype = 0; itype < result_dims[0]; ++itype) {
      for (npy_intp isize = 0; isize < result_dims[1]; ++isize) {
        for (npy_intp ilambda = 0; ilambda < result_dims[2]; ++ilambda) {
          for (npy_intp itheta = 0; itheta < result_dims[3]; ++itheta) {
            npy_intp iresult = 0;
            if (do_absorption) {
              // we use do_extinction as the result index, because the
              // absorption coefficients will either be result 0 or 1 depending
              // on whether extinction coefficients are present
              const npy_intp result_index = result_key->get_result_index(
                  itype, isize, ilambda, do_extinction);
              const AbsorptionCoefficientResult &result =
                  *static_cast<AbsorptionCoefficientResult *>(
                      results[result_index]);
              npy_intp Qabs_array_index[5] = {itype, isize, ilambda, itheta,
                                              iresult};
              ++iresult;
              npy_intp Qabspol_array_index[5] = {itype, isize, ilambda, itheta,
                                                 iresult};
              ++iresult;
              *reinterpret_cast<double *>(
                  PyArray_GetPtr(result_array, Qabs_array_index)) =
                  static_cast<double>(result.get_Qabs(itheta));
              *reinterpret_cast<double *>(
                  PyArray_GetPtr(result_array, Qabspol_array_index)) =
                  static_cast<double>(result.get_Qabspol(itheta));
            }
            if (do_extinction) {
              const npy_intp result_index =
                  result_key->get_result_index(itype, isize, ilambda, 0);
              const ExtinctionCoefficientResult &result =
                  *static_cast<ExtinctionCoefficientResult *>(
                      results[result_index]);
              npy_intp Qext_array_index[5] = {itype, isize, ilambda, itheta,
                                              iresult};
              ++iresult;
              npy_intp Qextpol_array_index[5] = {itype, isize, ilambda, itheta,
                                                 iresult};
              ++iresult;
              npy_intp Qextcpol_array_index[5] = {itype, isize, ilambda, itheta,
                                                  iresult};
              ++iresult;
              *reinterpret_cast<double *>(
                  PyArray_GetPtr(result_array, Qext_array_index)) =
                  static_cast<double>(result.get_Qext(itheta));
              *reinterpret_cast<double *>(
                  PyArray_GetPtr(result_array, Qextpol_array_index)) =
                  static_cast<double>(result.get_Qextpol(itheta));
              *reinterpret_cast<double *>(
                  PyArray_GetPtr(result_array, Qextcpol_array_index)) =
                  static_cast<double>(result.get_Qextcpol(itheta));
            }
          } // theta
        }   // lambda
      }     // size
    }       // type
    first_return_value = PyArray_Squeeze(result_array);
  }

  PyObject *second_return_value = nullptr;
  if (do_scattering) {
    npy_intp result_dims[7] = {
        static_cast<npy_intp>(result_key->composition_size()),
        static_cast<npy_intp>(result_key->size_size()),
        static_cast<npy_intp>(result_key->wavelength_size()),
        static_cast<npy_intp>(thetas.size()),
        static_cast<npy_intp>(2 * thetas.size()),
        4,
        4};
    PyArrayObject *result_array = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(7, result_dims, NPY_DOUBLE));

    uint_fast32_t num_results = do_extinction;
    num_results += do_absorption;
    for (npy_intp itype = 0; itype < result_dims[0]; ++itype) {
      for (npy_intp isize = 0; isize < result_dims[1]; ++isize) {
        for (npy_intp ilambda = 0; ilambda < result_dims[2]; ++ilambda) {
          const npy_intp result_index =
              result_key->get_result_index(itype, isize, ilambda, num_results);
          const ScatteringMatrixResult &result =
              *static_cast<ScatteringMatrixResult *>(results[result_index]);
          for (npy_intp itheta = 0; itheta < result_dims[3]; ++itheta) {
            for (npy_intp iphi = 0; iphi < result_dims[4]; ++iphi) {
              const Matrix<float_type> &Z = result.get_scattering_matrix(
                  2 * itheta * thetas.size() + iphi);
              for (npy_intp irow = 0; irow < 4; ++irow) {
                for (npy_intp icol = 0; icol < 4; ++icol) {
                  npy_intp array_index[7] = {itype, isize, ilambda, itheta,
                                             iphi,  irow,  icol};
                  *reinterpret_cast<double *>(
                      PyArray_GetPtr(result_array, array_index)) =
                      static_cast<double>(Z(irow, icol));
                } // col
              }   // row
            }     // phi
          }       // theta
        }         // lambda
      }           // size
    }             // type
    second_return_value = PyArray_Squeeze(result_array);
  }

  // now we are also done with the results, delete them too
  delete result_key;
  for (uint_fast32_t i = 0; i < results.size(); ++i) {
    delete results[i];
  }

  if (first_return_value != nullptr && second_return_value != nullptr) {
    return Py_BuildValue("OO", first_return_value, second_return_value);
  } else {
    if (first_return_value) {
      return first_return_value;
    } else {
      return second_return_value;
    }
  }
}

/*! @brief Methods exposed by the CTMmodule. */
static PyMethodDef CTMmethods[] = {
    {"get_equal_volume_radius",
     reinterpret_cast<PyCFunction>(get_equal_volume_radius),
     METH_VARARGS | METH_KEYWORDS,
     "Get the equal volume radius for the particle with the given equal area "
     "radius and axis ratio."},
    {"get_table", reinterpret_cast<PyCFunction>(get_table),
     METH_VARARGS | METH_KEYWORDS,
     "Test to see if the task-based algorithm can be coupled to the Python "
     "module."},
    {nullptr, nullptr, 0, nullptr}};

/*! @brief Module definition for the CTM module. */
static struct PyModuleDef CTMmodule = {PyModuleDef_HEAD_INIT, "CosTuuM",
                                       "T-matrix module.", -1, CTMmethods};

/**
 * @brief CosTuuM initialisation function.
 *
 * @return Pointer to the initialised module object.
 */
PyMODINIT_FUNC PyInit_CosTuuM() {

  // we need to call this to enable NumPy array functionality
  import_array();

  // create the module object
  PyObject *m = PyModule_Create(&CTMmodule);

  // old T-matrix object
  PyTMatrix::initialize(m);

  // add constants for the different orientation distribution types that can
  // be used for aligned grains
  PyModule_AddIntConstant(m, "DAVIS_GREENSTEIN_ALIGNMENT", 0);
  PyModule_AddIntConstant(m, "MISHCHENKO_ALIGNMENT", 1);
  PyModule_AddIntConstant(m, "DISABLE_ALIGNMENT", 2);

  // add constants for the different material types that can be used
  PyModule_AddIntConstant(m, "CARBON", 0);
  PyModule_AddIntConstant(m, "SILICON", 1);

  PySingleShapeShapeDistribution::initialize(m);
  PyDraineHensleyShapeDistribution::initialize(m);

  PySizeBasedAlignmentDistribution::initialize(m);

  PyDraineDustProperties::initialize(m);
  PyCustomDustProperties::initialize(m);

  PyOrientationDistribution::initialize(m);
  PyMishchenkoOrientationDistribution::initialize(m);
  PyDavisGreensteinOrientationDistribution::initialize(m);

  // return the module object
  return m;
}
