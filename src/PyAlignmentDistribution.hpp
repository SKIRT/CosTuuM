/**
 * @file PyAlignmentDistribution.hpp
 *
 * @brief Python exposure of the AlignmentDistribution interface and its
 * implementations.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYALIGNMENTDISTRIBUTION_HPP
#define PYALIGNMENTDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "SizeBasedAlignmentDistribution.hpp"

#include <Python.h>
#include <sstream>

/*! @brief Methods exposed to Python. */
static PyMethodDef AlignmentDistributionMethods[] = {{nullptr}};

/*! @brief Type for SizeBasedAlignmentDistribution wrapper objects. */
static PyTypeObject PySizeBasedAlignmentDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/**
 * @brief C struct wrapper around an AlignmentDistribution object.
 */
class PyAlignmentDistribution {
public:
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief Wrapped AlignmentDistribution object. */
  AlignmentDistribution *_alignment_distribution;

  /**
   * @brief Destructor for the AlignmentDistribution wrapper object.
   *
   * @param self PyAlignmentDistribution that is being deallocated.
   */
  static void dealloc(PyAlignmentDistribution *self) {
    // first free the wrapped AlignmentDistribution object
    // note that it doesn't matter which object we wrap here
    delete self->_alignment_distribution;
    // then free the Python object
    Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
  }

  /**
   * @brief Constructor for the AlignmentDistribution wrapper object.
   *
   * We use a single constructor for all AlignmentDistribution implementations
   * and initialise the wrapped object to a nullptr. The wrapped object is
   * created later, in the initialisation function, depending on the specific
   * implementation that is chosen.
   *
   * @param type Type for the object.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Pointer to the newly created object.
   */
  static PyObject *alloc(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    // allocate the Python object
    PyAlignmentDistribution *self =
        reinterpret_cast<PyAlignmentDistribution *>(type->tp_alloc(type, 0));
    // set the wrapped object to a nullptr, the object will be initialised
    // later, depending on the specific implementation that is chosen
    self->_alignment_distribution = nullptr;
    return reinterpret_cast<PyObject *>(self);
  }
};

/**
 * @brief C wrapper for the SizeBasedAlignmentDistribution implementation.
 */
class PySizeBasedAlignmentDistribution : public PyAlignmentDistribution {
public:
  /**
   * @brief init() function for a SizeBasedAlignmentDistribution object.
   *
   * Required arguments are:
   *  - minimum_size: Transition size above which dust grains are assumed to
   * align with the magnetic field (in m).
   *  - aligned_orientation_distribution_type: Type of orientation distribution
   * to use for aligned grains (needs proper documentation).
   *
   * Optional arguments are:
   *  - maximum_order: Maximum order of spherical basis function expansions
   *    (default: 100).
   *
   * @param self AlignmentDistribution wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyAlignmentDistribution *self, PyObject *args,
                  PyObject *kwargs) {

    // required arguments
    float_type minimum_size;
    int_fast32_t aligned_orientation_distribution_type;

    // optional arguments
    uint_fast32_t maximum_order = 100;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("minimum_size"),
                             strdup("aligned_orientation_distribution_type"),
                             strdup("maximum_order"), nullptr};

    // allocate temporary variables to store double precision arguments
    double minimum_size_d;
    // temporary variables for integer arguments
    int aligned_orientation_distribution_type_i;
    unsigned int maximum_order_i = maximum_order;
    // parse the keywords/positional arguments
    // d is a double
    // I is an integer
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "di|I", kwlist, &minimum_size_d,
            &aligned_orientation_distribution_type_i, &maximum_order_i)) {
      // PyArg_ParseTupleAndKeywords will return 0 if a required argument was
      // missing, if an argument of the wrong type was provided or if the number
      // of arguments does not match the expectation
      // we use ctm_warning so that we can gracefully exit and let Python handle
      // the exception
      // ctm_error would call abort, which would kill the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // exit code 1 signals to Python that something was wrong
      return 1;
    }
    // unpack double precision arguments
    minimum_size = minimum_size_d;
    // unpack integer arguments
    aligned_orientation_distribution_type =
        aligned_orientation_distribution_type_i;
    maximum_order = maximum_order_i;

    // create the object
    self->_alignment_distribution = new SizeBasedAlignmentDistribution(
        minimum_size, aligned_orientation_distribution_type, maximum_order);

    return 0;
  }

  /**
   * @brief Initialize the PySizeBasedAlignmentDistribution object and add it to
   * the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".SizeBasedAlignmentDistribution";
    // set the required fields in the ObjectType struct
    PySizeBasedAlignmentDistributionType.tp_name = object_name.str().c_str();
    PySizeBasedAlignmentDistributionType.tp_basicsize =
        sizeof(PySizeBasedAlignmentDistribution);
    PySizeBasedAlignmentDistributionType.tp_dealloc = (destructor)dealloc;
    PySizeBasedAlignmentDistributionType.tp_methods =
        AlignmentDistributionMethods;
    PySizeBasedAlignmentDistributionType.tp_init = (initproc)init;
    PySizeBasedAlignmentDistributionType.tp_new = alloc;

    // finalize creation of the PySizeBasedAlignmentDistributionType
    PyType_Ready(&PySizeBasedAlignmentDistributionType);
    // add a PySizeBasedAlignmentDistributionType to the module
    Py_INCREF(&PySizeBasedAlignmentDistributionType);
    PyModule_AddObject(
        module, "SizeBasedAlignmentDistribution",
        reinterpret_cast<PyObject *>(&PySizeBasedAlignmentDistributionType));
  }
};

#endif // PYALIGNMENTDISTRIBUTION_HPP
