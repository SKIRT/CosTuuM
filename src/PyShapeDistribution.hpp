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
 * @file PyShapeDistribution.hpp
 *
 * @brief Python exposure of the ShapeDistribution interface and its
 * implementations.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYSHAPEDISTRIBUTION_HPP
#define PYSHAPEDISTRIBUTION_HPP

#include "Configuration.hpp"
#include "DraineHensleyShapeDistribution.hpp"
#include "SingleShapeShapeDistribution.hpp"

#include <Python.h>
#include <sstream>

/// ShapeDistribution object

/*! @brief Type for SingleShapeShapeDistribution wrapper objects. */
static PyTypeObject PySingleShapeShapeDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/*! @brief Type for DraineHensleyShapeDistribution wrapper objects. */
static PyTypeObject PyDraineHensleyShapeDistributionType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/**
 * @brief C struct wrapper around a ShapeDistribution object.
 */
class PyShapeDistribution {
public:
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief Wrapped ShapeDistribution object. */
  ShapeDistribution *_shape_distribution;

  /**
   * @brief Destructor for the ShapeDistribution wrapper object.
   *
   * @param self ShapeDistributionObject that is being deallocated.
   */
  static void dealloc(PyShapeDistribution *self) {
    // first free the wrapped ShapeDistribution object
    // note that it doesn't matter which object we wrap here
    delete self->_shape_distribution;
    // then free the Python object
    Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
  }

  /**
   * @brief Constructor for the ShapeDistribution wrapper object.
   *
   * We use a single constructor for all ShapeDistribution implementations and
   * initialise the wrapped object to a nullptr. The wrapped object is created
   * later, in the initialisation function, depending on the specific
   * implementation that is chosen.
   *
   * @param type Type for the object.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Pointer to the newly created object.
   */
  static PyObject *alloc(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    // allocate the Python object
    PyShapeDistribution *self =
        reinterpret_cast<PyShapeDistribution *>(type->tp_alloc(type, 0));
    // set the wrapped object to a nullptr, the object will be initialised
    // later, depending on the specific implementation that is chosen
    self->_shape_distribution = nullptr;
    return reinterpret_cast<PyObject *>(self);
  }

  /**
   * @brief Get the limits for the shape distribution.
   *
   * @param self PyDustProperties object.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Nothing.
   */
  static PyObject *get_limits(PyShapeDistribution *self, PyObject *args,
                              PyObject *kwargs) {

    return Py_BuildValue(
        "dd",
        static_cast<double>(
            self->_shape_distribution->get_minimum_axis_ratio()),
        static_cast<double>(
            self->_shape_distribution->get_maximum_axis_ratio()));
  }
};

/*! @brief Methods exposed to Python. */
static PyMethodDef ShapeDistributionMethods[] = {
    {"get_limits",
     reinterpret_cast<PyCFunction>(PyShapeDistribution::get_limits),
     METH_VARARGS | METH_KEYWORDS,
     "Get the limits for the shape distribution."},
    {nullptr}};

/**
 * @brief C wrapper for the SingleShapeShapeDistribution implementation.
 */
class PySingleShapeShapeDistribution : public PyShapeDistribution {
public:
  /**
   * @brief init() function for a SingleShapeShapeDistribution object.
   *
   * Required arguments are:
   *  - axis_ratio: Axis ratio parameter for the single shape contained in this
   *    distribution.
   *
   * @param self ShapeDistribution wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyShapeDistribution *self, PyObject *args, PyObject *kwargs) {

    // required arguments
    float_type axis_ratio;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("axis_ratio"), nullptr};

    // allocate temporary variables to store double precision arguments
    double axis_ratio_d;
    // parse the keywords/positional arguments
    // d is a double
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist,
                                     &axis_ratio_d)) {
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

    // create the object
    self->_shape_distribution = new SingleShapeShapeDistribution(axis_ratio);

    return 0;
  }

  /**
   * @brief Initialize the PySingleShapeShapeDistribution object and add it to
   * the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".SingleShapeShapeDistribution";
    // set the required fields in the ObjectType struct
    PySingleShapeShapeDistributionType.tp_name = object_name.str().c_str();
    PySingleShapeShapeDistributionType.tp_basicsize =
        sizeof(PySingleShapeShapeDistribution);
    PySingleShapeShapeDistributionType.tp_dealloc = (destructor)dealloc;
    PySingleShapeShapeDistributionType.tp_methods = ShapeDistributionMethods;
    PySingleShapeShapeDistributionType.tp_init = (initproc)init;
    PySingleShapeShapeDistributionType.tp_new = alloc;

    // finalize creation of the PySingleShapeShapeDistributionType
    PyType_Ready(&PySingleShapeShapeDistributionType);
    // add a PySingleShapeShapeDistributionType to the module
    Py_INCREF(&PySingleShapeShapeDistributionType);
    PyModule_AddObject(
        module, "SingleShapeShapeDistribution",
        reinterpret_cast<PyObject *>(&PySingleShapeShapeDistributionType));
  }
};

/**
 * @brief C wrapper for the DraineHensleyShapeDistribution implementation.
 */
class PyDraineHensleyShapeDistribution : public PyShapeDistribution {
public:
  /**
   * @brief init() function for a DraineHensleyShapeDistribution object.
   *
   * Optional arguments are:
   *  - numpoints: Number of points to use to sample the distribution (default:
   *    20).
   *  - cutoff: Probability value that determines the lower and upper bounds
   *    for the sampling interval of the distribution (default: 0.15).
   *  - fraction: Fraction of the distribution that is sampled, both the oblate
   *    and prolate shapes are sampled with the same fraction. This parameter
   *    is incompatible with the cutoff parameter and is disabled by default
   *    (default: -1).
   *
   * @param self ShapeDistribution wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyShapeDistribution *self, PyObject *args, PyObject *kwargs) {

    // optional arguments
    uint_fast32_t npoints = 20;
    float_type cutoff = 0.15;
    float_type fraction = -1.;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("npoints"), strdup("cutoff"),
                             strdup("fraction"), nullptr};

    // placeholders for integer arguments
    unsigned int npoints_i = npoints;
    // placeholders for float arguments
    double cutoff_d = static_cast<double>(cutoff);
    double fraction_d = static_cast<double>(fraction);
    // parse the keywords/positional arguments
    // I is an unsigned integer
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|Idd", kwlist, &npoints_i,
                                     &cutoff_d, &fraction_d)) {
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
    npoints = npoints_i;
    // unpack float arguments
    cutoff = cutoff_d;
    fraction = fraction_d;

    // create the object
    self->_shape_distribution =
        new DraineHensleyShapeDistribution(npoints, cutoff, fraction);

    return 0;
  }

  /**
   * @brief Initialize the PyDraineHensleyShapeDistribution object and add it to
   * the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".DraineHensleyShapeDistribution";
    // set the required fields in the ObjectType struct
    PyDraineHensleyShapeDistributionType.tp_name = object_name.str().c_str();
    PyDraineHensleyShapeDistributionType.tp_basicsize =
        sizeof(PyDraineHensleyShapeDistribution);
    PyDraineHensleyShapeDistributionType.tp_dealloc = (destructor)dealloc;
    PyDraineHensleyShapeDistributionType.tp_methods = ShapeDistributionMethods;
    PyDraineHensleyShapeDistributionType.tp_init = (initproc)init;
    PyDraineHensleyShapeDistributionType.tp_new = alloc;

    // finalize creation of the PySingleShapeShapeDistributionType
    PyType_Ready(&PyDraineHensleyShapeDistributionType);
    // add a PySingleShapeShapeDistributionType to the module
    Py_INCREF(&PyDraineHensleyShapeDistributionType);
    PyModule_AddObject(
        module, "DraineHensleyShapeDistribution",
        reinterpret_cast<PyObject *>(&PyDraineHensleyShapeDistributionType));
  }
};

#endif // PYSHAPEDISTRIBUTION_HPP
