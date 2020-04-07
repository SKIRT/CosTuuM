/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file PyOrientationDistribution.hpp
 *
 * @brief Python exposure of the OrientationDistribution interface and its
 * implementations.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYORIENTATIONDISTRIBUTION_HPP
#define PYORIENTATIONDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "PyHelperFunctions.hpp"

#include <Python.h>
#include <sstream>

/*! @brief Type for normal OrientationDistribution wrapper objects. */
static PyTypeObject RandomOrientationDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/*! @brief Type for MishchenkoOrientationDistribution wrapper objects. */
static PyTypeObject MishchenkoOrientationDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

///*! @brief Type for DavisGreensteinOrientationDistribution wrapper objects. */
static PyTypeObject DavisGreensteinOrientationDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

///*! @brief Type for CustomOrientationDistribution wrapper objects. */
static PyTypeObject CustomOrientationDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/**
 * @brief C struct wrapper around an OrientationDistribution object.
 */
class PyOrientationDistribution {
public:
  /*! @brief Python object members. */
  PyObject_HEAD;

  /*! @brief Wrapped OrientationDistribution object. */
  OrientationDistribution *_orientation_distribution;

  /**
   * @brief Destructor for the OrientationDistribution wrapper object.
   *
   * @param self PyOrientationDistribution that is being deallocated.
   */
  static void dealloc(PyOrientationDistribution *self) {
    // first free the wrapped AlignmentDistribution object
    // note that it doesn't matter which object we wrap here
    delete self->_orientation_distribution;
    // then free the Python object
    Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
  }

  /**
   * @brief Constructor for the OrientationDistribution wrapper object.
   *
   * We use a single constructor for all OrientationDistribution implementations
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
    PyOrientationDistribution *self =
        reinterpret_cast<PyOrientationDistribution *>(type->tp_alloc(type, 0));
    // set the wrapped object to a nullptr, the object will be initialised
    // later, depending on the specific implementation that is chosen
    self->_orientation_distribution = nullptr;
    return reinterpret_cast<PyObject *>(self);
  }

  /**
   * @brief init() function for a RandomOrientationDistribution object.
   *
   * Optional arguments are:
   *  - order: Order of coefficient expansion for the distribution
   *    (default: 100).
   *
   * @param self OrientationDistribution wrapper object that is being
   * initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyOrientationDistribution *self, PyObject *args,
                  PyObject *kwargs) {

    // optional arguments
    uint_fast32_t order = 100;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("order"), nullptr};

    // temporary variables for integer arguments
    unsigned int order_i = order;
    // parse the keywords/positional arguments
    // I is an integer
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|I", kwlist, &order_i)) {
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
    // unpack integer arguments
    order = order_i;

    // create the object
    self->_orientation_distribution = new OrientationDistribution(order);
    self->_orientation_distribution->initialise();

    return 0;
  }

  /**
   * @brief Get the value of the orientation distribution for the given input
   * zenith angle(s).
   *
   * Required arguments are:
   *  - theta: Input zenith angle(s) (in radians).
   *
   * @param self PyOrientationDistribution object.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Nothing.
   */
  static PyObject *evaluate(PyOrientationDistribution *self, PyObject *args,
                            PyObject *kwargs) {

    // parse arguments //

    // required arguments
    PyArrayObject *input_thetas;

    // list of keywords
    static char *kwlist[] = {strdup("theta"), nullptr};

    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&", kwlist,
                                     PyArray_Converter, &input_thetas)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }

    std::vector<float_type> thetas =
        unpack_numpy_array<float_type, double>(input_thetas);
    Py_DECREF(input_thetas);

    npy_intp result_dims[1] = {static_cast<npy_intp>(thetas.size())};
    PyArrayObject *result_array = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(1, result_dims, NPY_DOUBLE));

    const OrientationDistribution &distribution =
        *self->_orientation_distribution;
    for (npy_intp i = 0; i < result_dims[0]; ++i) {

      std::vector<float_type> wigner_d(distribution.get_maximum_order(), 0.);
      std::vector<float_type> dwigner_d(distribution.get_maximum_order(), 0.);
      // we have to manually compute the zeroth order value
      wigner_d[0] = 1.;
      SpecialFunctions::wigner_dn_0m(cos(thetas[i]),
                                     distribution.get_maximum_order() - 1, 0,
                                     &wigner_d[1], &dwigner_d[1]);

      float_type value = 0.;
      for (uint_fast32_t n = 0; n < distribution.get_maximum_order(); ++n) {
        value +=
            0.5 * (2. * n + 1.) * distribution.get_coefficient(n) * wigner_d[n];
      }

      *reinterpret_cast<double *>(PyArray_GETPTR1(result_array, i)) =
          static_cast<double>(value);
    }

    return PyArray_Squeeze(result_array);
  }

  /**
   * @brief Initialize the RandomOrientationDistributionType object and add it
   * to the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module);
};

/*! @brief Methods exposed to Python. */
static PyMethodDef OrientationDistributionMethods[] = {
    {"evaluate",
     reinterpret_cast<PyCFunction>(PyOrientationDistribution::evaluate),
     METH_VARARGS | METH_KEYWORDS,
     "Get the value of the orientation distribution for the given input zenith "
     "angle(s)."},
    {nullptr}};

void PyOrientationDistribution::initialize(PyObject *module) {

  std::stringstream object_name;
  object_name << PyModule_GetName(module);
  object_name << ".RandomOrientationDistribution";
  // set the required fields in the ObjectType struct
  RandomOrientationDistributionType.tp_name = object_name.str().c_str();
  RandomOrientationDistributionType.tp_basicsize =
      sizeof(PyOrientationDistribution);
  RandomOrientationDistributionType.tp_dealloc = (destructor)dealloc;
  RandomOrientationDistributionType.tp_methods = OrientationDistributionMethods;
  RandomOrientationDistributionType.tp_init = (initproc)init;
  RandomOrientationDistributionType.tp_new = alloc;

  // finalize creation of the RandomOrientationDistributionType
  PyType_Ready(&RandomOrientationDistributionType);
  // add a RandomOrientationDistributionType to the module
  Py_INCREF(&RandomOrientationDistributionType);
  PyModule_AddObject(
      module, "RandomOrientationDistribution",
      reinterpret_cast<PyObject *>(&RandomOrientationDistributionType));
}

/**
 * @brief C wrapper for the MishchenkoOrientationDistribution implementation.
 */
class PyMishchenkoOrientationDistribution : public PyOrientationDistribution {
public:
  /**
   * @brief init() function for a MishchenkoOrientationDistribution object.
   *
   * Required arguments are:
   *  - cos2beta: Alignment strength parameter.
   *
   * Optional arguments are:
   *  - order: Order of coefficient expansion for the distribution
   *    (default: 100).
   *
   * @param self MishchenkoOrientationDistribution wrapper object that is being
   * initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyMishchenkoOrientationDistribution *self, PyObject *args,
                  PyObject *kwargs) {

    // required arguments
    float_type cos2beta;

    // optional arguments
    uint_fast32_t order = 100;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("cos2beta"), strdup("order"), nullptr};

    // allocate temporary variables to store double precision arguments
    double cos2beta_d;
    // temporary variables for integer arguments
    unsigned int order_i = order;
    // parse the keywords/positional arguments
    // d is a double
    // I is an integer
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|I", kwlist, &cos2beta_d,
                                     &order_i)) {
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
    cos2beta = cos2beta_d;
    // unpack integer arguments
    order = order_i;

    // create the object
    self->_orientation_distribution =
        new MishchenkoOrientationDistribution(order, cos2beta);

    return 0;
  }

  /**
   * @brief Initialize the PyMishchenkoOrientationDistribution object and add it
   * to the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".MishchenkoOrientationDistribution";
    // set the required fields in the ObjectType struct
    MishchenkoOrientationDistributionType.tp_name = object_name.str().c_str();
    MishchenkoOrientationDistributionType.tp_basicsize =
        sizeof(PyMishchenkoOrientationDistribution);
    MishchenkoOrientationDistributionType.tp_dealloc = (destructor)dealloc;
    MishchenkoOrientationDistributionType.tp_methods =
        OrientationDistributionMethods;
    MishchenkoOrientationDistributionType.tp_init = (initproc)init;
    MishchenkoOrientationDistributionType.tp_new = alloc;

    // finalize creation of the MishchenkoOrientationDistributionType
    PyType_Ready(&MishchenkoOrientationDistributionType);
    // add a MishchenkoOrientationDistributionType to the module
    Py_INCREF(&MishchenkoOrientationDistributionType);
    PyModule_AddObject(
        module, "MishchenkoOrientationDistribution",
        reinterpret_cast<PyObject *>(&MishchenkoOrientationDistributionType));
  }
};

/**
 * @brief C wrapper for the DavisGreensteinOrientationDistribution
 * implementation.
 */
class PyDavisGreensteinOrientationDistribution
    : public PyOrientationDistribution {
public:
  /**
   * @brief init() function for a DavisGreensteinOrientationDistribution object.
   *
   * Required arguments are:
   *  - axis_ratio: Axis ratio of the particles that align (to distinguish
   *    between oblate and prolate).
   *
   * Optional arguments are:
   *  - order: Order of coefficient expansion for the distribution
   *    (default: 100).
   *
   * @param self DavisGreensteinOrientationDistribution wrapper object that is
   * being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyDavisGreensteinOrientationDistribution *self,
                  PyObject *args, PyObject *kwargs) {

    // required arguments
    float_type axis_ratio;

    // optional arguments
    uint_fast32_t order = 100;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("axis_ratio"), strdup("order"), nullptr};

    // allocate temporary variables to store double precision arguments
    double axis_ratio_d;
    // temporary variables for integer arguments
    unsigned int order_i = order;
    // parse the keywords/positional arguments
    // d is a double
    // I is an integer
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|I", kwlist, &axis_ratio_d,
                                     &order_i)) {
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
    axis_ratio = axis_ratio_d;
    // unpack integer arguments
    order = order_i;

    // create the object
    self->_orientation_distribution =
        new DavisGreensteinOrientationDistribution(order, axis_ratio);

    return 0;
  }

  /**
   * @brief Initialize the PyDavisGreensteinOrientationDistribution object and
   * add it to the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".DavisGreensteinOrientationDistribution";
    // set the required fields in the ObjectType struct
    DavisGreensteinOrientationDistributionType.tp_name =
        object_name.str().c_str();
    DavisGreensteinOrientationDistributionType.tp_basicsize =
        sizeof(PyDavisGreensteinOrientationDistribution);
    DavisGreensteinOrientationDistributionType.tp_dealloc = (destructor)dealloc;
    DavisGreensteinOrientationDistributionType.tp_methods =
        OrientationDistributionMethods;
    DavisGreensteinOrientationDistributionType.tp_init = (initproc)init;
    DavisGreensteinOrientationDistributionType.tp_new = alloc;

    // finalize creation of the DavisGreensteinOrientationDistributionType
    PyType_Ready(&DavisGreensteinOrientationDistributionType);
    // add a DavisGreensteinOrientationDistributionType to the module
    Py_INCREF(&DavisGreensteinOrientationDistributionType);
    PyModule_AddObject(module, "DavisGreensteinOrientationDistribution",
                       reinterpret_cast<PyObject *>(
                           &DavisGreensteinOrientationDistributionType));
  }
};

/**
 * @brief OrientationDistribution implementation that uses a user-provided
 * distribution function.
 */
class CustomOrientationDistribution : public OrientationDistribution {
private:
  /*! @brief Python distribution function. */
  PyObject *_python_function;

  /*! @brief Additional arguments for the Python distribution function. */
  PyObject *_additional_arguments;

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Highest order for which we store a coefficient.
   * @param python_function Python distribution function.
   * @param additional_arguments Additional arguments for the Python function.
   */
  inline CustomOrientationDistribution(const uint_fast32_t nmax,
                                       PyObject *python_function,
                                       PyObject *additional_arguments = nullptr)
      : OrientationDistribution(nmax), _python_function(python_function),
        _additional_arguments(additional_arguments) {
    // tell Python we store a copy of the function object
    Py_INCREF(_python_function);
    if (_additional_arguments) {
      Py_INCREF(_additional_arguments);
    }
    initialise(1.e-4, 1.e-4);
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~CustomOrientationDistribution() {
    // tell Python we are done with our local copy of the function object
    Py_DECREF(_python_function);
    if (_additional_arguments) {
      Py_DECREF(_additional_arguments);
    }
  }

  /**
   * @brief Virtual function that contains an actual expression for the
   * orientation distribution function as a function of the zenith angle
   * @f$\beta{}@f$.
   *
   * The zenith angle @f$\beta{}@f$ represents the angle between the average
   * orientation of the rotation axis of the spheroidal dust particles and the
   * specific orientation for which we want to know the probability.
   * A distribution function @f$p(\beta{}) = \delta{}(\beta{})@f$ e.g. would
   * correspond to perfect alignment of all dust particles.
   *
   * The additional cosine and sine arguments are provided because they often
   * appear in more complicated expressions for @f$p(\beta{})@f$, and also
   * because they are computed in the integrand anyway.
   *
   * @param beta Zenith angle @f$\beta{} \in{} [0, \pi{}]@f$.
   * @param cosbeta Cosine of the zenith angle, @f$\cos(\beta{})@f$.
   * @param sinbeta Sine of the zenith angle, @f$\sin(\beta{})@f$.
   * @return Orientation distribution function for that value of
   * @f$\beta{}@f$, @f$p(\beta{})@f$.
   */
  virtual float_type operator()(const float_type beta, const float_type cosbeta,
                                const float_type sinbeta) const {

    // wrappers for float_type values
    const double beta_d = static_cast<double>(beta);
    // build the argument list for the Python function call
    PyObject *arglist;

    if (_additional_arguments) {
      arglist = Py_BuildValue("(dO)", beta_d, _additional_arguments);
    } else {
      arglist = Py_BuildValue("(d)", beta_d);
    }
    // place the Python function call
    PyObject *result = PyEval_CallObject(_python_function, arglist);
    // we are done with the argument list
    Py_DECREF(arglist);

    if (result == nullptr) {
      ctm_warning("Error during function call!");
    }

    // parse the result
    double result_value;
    // first check what kind of result we got
    if (!PyFloat_Check(result)) {
      // we did not get a float. Check if we got a NumPy array instead
      if (PyArray_Check(result)) {
        // we did
        PyArrayObject *array = reinterpret_cast<PyArrayObject *>(result);
        // check if we got a scalar array
        if (!PyArray_CheckScalar(array)) {
          // nope. Bail out.
          ctm_error("Function returns an array!");
        }
        // check if the array has float elements
        if (!PyArray_ISFLOAT(array)) {
          // nope. That doesn't work.
          ctm_error("Function does not return a float!");
        }
        // get the value by casting the array to a double
        double *results = reinterpret_cast<double *>(PyArray_BYTES(array));
        result_value = results[0];
      } else {
        // nope. Give up.
        ctm_error("Function does not return a float!");
      }
    } else {
      // get the real and imaginary parts of the complex number
      result_value = PyFloat_AsDouble(result);
    }

    // we are done with the result
    Py_DECREF(result);

    return result_value;
  }
};

/**
 * @brief C wrapper for the CustomOrientationDistribution implementation.
 */
class PyCustomOrientationDistribution : public PyOrientationDistribution {
public:
  /**
   * @brief init() function for a CustomOrientationDistribution object.
   *
   * Required arguments are:
   *  - function: Python function that returns the value of the orientation
   *    distribution for a given input zenith angle.
   *
   * Optional arguments are:
   *  - order: Order of coefficient expansion for the distribution
   *    (default: 100).
   *  - additional_arguments: Additional arguments for the Python function.
   *
   * @param self OrientationDistribution wrapper object that is being
   * initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyOrientationDistribution *self, PyObject *args,
                  PyObject *kwargs) {

    // required arguments
    PyObject *function;

    // optional arguments
    uint_fast32_t order = 100;
    PyObject *additional_arguments = nullptr;

    // list of keywords
    static char *kwlist[] = {strdup("function"), strdup("order"),
                             strdup("additional_arguments"), nullptr};

    // temporary variables for integer arguments
    unsigned int order_i = order;
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|IO", kwlist, &function,
                                     &order_i, &additional_arguments)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      return 1;
    }
    // unpack integer arguments
    order = order_i;

    self->_orientation_distribution = new CustomOrientationDistribution(
        order, function, additional_arguments);

    return 0;
  }

  /**
   * @brief Initialize the PyCustomOrientationDistribution object and add it to
   * the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".CustomOrientationDistribution";
    // set the required fields in the ObjectType struct
    CustomOrientationDistributionType.tp_name = object_name.str().c_str();
    CustomOrientationDistributionType.tp_basicsize =
        sizeof(PyCustomOrientationDistribution);
    CustomOrientationDistributionType.tp_dealloc = (destructor)dealloc;
    CustomOrientationDistributionType.tp_methods =
        OrientationDistributionMethods;
    CustomOrientationDistributionType.tp_init = (initproc)init;
    CustomOrientationDistributionType.tp_new = alloc;

    // finalize creation of the CustomOrientationDistributionType
    PyType_Ready(&CustomOrientationDistributionType);
    // add a CustomOrientationDistributionType to the module
    Py_INCREF(&CustomOrientationDistributionType);
    PyModule_AddObject(
        module, "CustomOrientationDistribution",
        reinterpret_cast<PyObject *>(&CustomOrientationDistributionType));
  }
};

#endif // PYORIENTATIONDISTRIBUTION_HPP
