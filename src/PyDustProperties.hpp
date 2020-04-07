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
 * @file PyDustProperties.hpp
 *
 * @brief Python exposure of the DustProperties interface and its
 * implementations.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYDUSTPROPERTIES_HPP
#define PYDUSTPROPERTIES_HPP

#include "Configuration.hpp"
#include "DraineDustProperties.hpp"

#include <Python.h>
#include <cstring>
#include <sstream>

/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

/*! @brief Type for DraineDustProperties wrapper objects. */
static PyTypeObject PyDraineDustPropertiesType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/*! @brief Type for CustomDustProperties wrapper objects. */
static PyTypeObject PyCustomDustPropertiesType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/**
 * @brief C struct wrapper around a DustProperties object.
 */
class PyDustProperties {
public:
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief Wrapped DustProperties object. */
  DustProperties *_dust_properties;

  /**
   * @brief Destructor for the DustProperties wrapper object.
   *
   * @param self PyDustProperties that is being deallocated.
   */
  static void dealloc(PyDustProperties *self) {
    // first free the wrapped DustProperties object
    // note that it doesn't matter which object we wrap here
    delete self->_dust_properties;
    // then free the Python object
    Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
  }

  /**
   * @brief Constructor for the DustProperties wrapper object.
   *
   * We use a single constructor for all DustProperties implementations
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
    PyDustProperties *self =
        reinterpret_cast<PyDustProperties *>(type->tp_alloc(type, 0));
    // set the wrapped object to a nullptr, the object will be initialised
    // later, depending on the specific implementation that is chosen
    self->_dust_properties = nullptr;
    return reinterpret_cast<PyObject *>(self);
  }

  /**
   * @brief Get the refractive index for the given type, size and wavelength.
   *
   * Required arguments are:
   *  - type: Dust grain type.
   *  - size: Dust grain size (in m).
   *  - wavelength: Incoming/outgoing photon wavelength (in m).
   *
   * @param self PyDustProperties object.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Nothing.
   */
  static PyObject *get_refractive_index(PyDustProperties *self, PyObject *args,
                                        PyObject *kwargs) {

    // parse arguments //

    // required arguments
    int_fast32_t type;
    float_type size;
    float_type wavelength;

    // list of keywords
    static char *kwlist[] = {strdup("type"), strdup("size"),
                             strdup("wavelength"), nullptr};

    // placeholders for float_type arguments
    double size_d, wavelength_d;
    // placeholders for int arguments
    int type_i;
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "idd", kwlist, &type_i,
                                     &size_d, &wavelength_d)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }
    // convert float_type arguments
    size = size_d;
    wavelength = wavelength_d;
    // convert int arguments
    type = type_i;

    const std::complex<float_type> refractive_index =
        self->_dust_properties->get_refractive_index(wavelength, size, type);

    return PyComplex_FromDoubles(static_cast<double>(refractive_index.real()),
                                 static_cast<double>(refractive_index.imag()));
  }
};

/*! @brief Methods exposed to Python. */
static PyMethodDef DustPropertiesMethods[] = {
    {"get_refractive_index",
     reinterpret_cast<PyCFunction>(PyDustProperties::get_refractive_index),
     METH_VARARGS | METH_KEYWORDS,
     "Get the refractive index for the dust grain with the given properties "
     "and for the given wavelength."},
    {nullptr}};

/**
 * @brief C wrapper for the DraineDustProperties implementation.
 */
class PyDraineDustProperties : public PyDustProperties {
public:
  /**
   * @brief init() function for a DraineDustProperties object.
   *
   * Optional arguments are:
   *  - dust_temperature: Temperature of the dust grains (in K; default: 20. K).
   *
   * @param self DustProperties wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyDustProperties *self, PyObject *args, PyObject *kwargs) {

    // optional arguments
    float_type dust_temperature = 20.;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("dust_temperature"), nullptr};

    // allocate temporary variables to store double precision arguments
    double dust_temperature_d = static_cast<double>(dust_temperature);
    // parse the keywords/positional arguments
    // d is a double
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|d", kwlist,
                                     &dust_temperature_d)) {
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
    dust_temperature = dust_temperature_d;

    self->_dust_properties = new DraineDustProperties(dust_temperature);

    return 0;
  }

  /**
   * @brief Initialize the PyDraineDustProperties object and add it to the given
   * module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".DraineDustProperties";
    // set the required fields in the ObjectType struct
    PyDraineDustPropertiesType.tp_name = object_name.str().c_str();
    PyDraineDustPropertiesType.tp_basicsize = sizeof(PyDraineDustProperties);
    PyDraineDustPropertiesType.tp_dealloc = (destructor)dealloc;
    PyDraineDustPropertiesType.tp_methods = DustPropertiesMethods;
    PyDraineDustPropertiesType.tp_init = (initproc)init;
    PyDraineDustPropertiesType.tp_new = alloc;

    // finalize creation of the PyDraineDustPropertiesType
    PyType_Ready(&PyDraineDustPropertiesType);
    // add a PyDraineDustPropertiesType to the module
    Py_INCREF(&PyDraineDustPropertiesType);
    PyModule_AddObject(
        module, "DraineDustProperties",
        reinterpret_cast<PyObject *>(&PyDraineDustPropertiesType));
  }
};

/**
 * @brief DustProperties implementation that reads values from a user-provided
 * Python function.
 */
class CustomDustProperties : public DustProperties {
private:
  /*! @brief Python function that is being wrapped. */
  PyObject *_python_function;

  /*! @brief Additional arguments for the Python function. */
  PyObject *_additional_arguments;

public:
  /**
   * @brief Constructor.
   *
   * @param python_function Python function that is being wrapped.
   * @param additional_arguments Additional arguments for the Python function.
   */
  inline CustomDustProperties(PyObject *python_function,
                              PyObject *additional_arguments = nullptr)
      : _python_function(python_function),
        _additional_arguments(additional_arguments) {
    // tell Python we store a copy of the function object
    Py_INCREF(_python_function);
    if (_additional_arguments) {
      Py_INCREF(_additional_arguments);
    }
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~CustomDustProperties() {
    // tell Python we are done with our local copy of the function object
    Py_DECREF(_python_function);
    if (_additional_arguments) {
      Py_DECREF(_additional_arguments);
    }
  }

  /**
   * @brief Get the refractive index for the given wavelength, grain size and
   * grain type.
   *
   * @param wavelength Wavelength of the incoming radiation (in m).
   * @param grain_size Size of the grain (in m).
   * @param grain_type Type of dust grain.
   * @return Complex refractive index of the grain for radiation at this
   * wavelength.
   */
  virtual std::complex<float_type>
  get_refractive_index(const float_type wavelength, const float_type grain_size,
                       const int_fast32_t grain_type) const {

    // wrappers for float_type values
    const double wavelength_d = static_cast<double>(wavelength);
    const double grain_size_d = static_cast<double>(grain_size);
    // wrappers for int values
    const int grain_type_i = grain_type;
    // build the argument list for the Python function call
    PyObject *arglist;

    if (_additional_arguments) {
      arglist = Py_BuildValue("(ddiO)", wavelength_d, grain_size_d,
                              grain_type_i, _additional_arguments);
    } else {
      arglist =
          Py_BuildValue("(ddi)", wavelength_d, grain_size_d, grain_type_i);
    }
    // place the Python function call
    PyObject *result = PyEval_CallObject(_python_function, arglist);
    // we are done with the argument list
    Py_DECREF(arglist);

    if (result == nullptr) {
      ctm_warning("Error during function call!");
    }

    // parse the result
    double result_real, result_imag;
    // first check what kind of result we got
    if (strcmp(result->ob_type->tp_name, "complex") != 0) {
      // we did not get a complex number. Check if we got a NumPy array instead
      if (PyArray_Check(result)) {
        // we did
        PyArrayObject *array = reinterpret_cast<PyArrayObject *>(result);
        // check if we got a scalar array
        if (!PyArray_CheckScalar(array)) {
          // nope. Bail out.
          ctm_error("Function returns an array!");
        }
        // check if the array has complex elements
        if (!PyArray_ISCOMPLEX(array)) {
          // nope. That doesn't work.
          ctm_error("Function does not return a complex number!");
        }
        // get the real and imaginary parts by casting the array to two doubles
        double *results = reinterpret_cast<double *>(PyArray_BYTES(array));
        result_real = results[0];
        result_imag = results[1];
      } else {
        // nope. Give up.
        ctm_error("Function does not return a complex number!");
      }
    } else {
      // get the real and imaginary parts of the complex number
      result_real = PyComplex_RealAsDouble(result);
      result_imag = PyComplex_ImagAsDouble(result);
    }

    // we are done with the result
    Py_DECREF(result);

    return std::complex<float_type>(result_real, result_imag);
  }
};

/**
 * @brief C wrapper for the CustomDustProperties implementation.
 */
class PyCustomDustProperties : public PyDustProperties {
public:
  /**
   * @brief init() function for a CustomDustProperties object.
   *
   * Required arguments are:
   *  - function: Python function that returns the complex refractive index for
   *  a given input wavelength, dust grain size and dust grain type (in this
   *  order).
   *
   * Optional arguments are:
   *  - additional_arguments: Additional arguments for the Python function.
   *
   * @param self DustProperties wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyDustProperties *self, PyObject *args, PyObject *kwargs) {

    // required arguments
    PyObject *function;

    // optional arguments
    PyObject *additional_arguments = nullptr;

    // list of keywords
    static char *kwlist[] = {strdup("function"), strdup("additional_arguments"),
                             nullptr};

    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist, &function,
                                     &additional_arguments)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      return 1;
    }

    self->_dust_properties =
        new CustomDustProperties(function, additional_arguments);

    return 0;
  }

  /**
   * @brief Initialize the PyCustomDustProperties object and add it to the given
   * module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module) {

    std::stringstream object_name;
    object_name << PyModule_GetName(module);
    object_name << ".CustomDustProperties";
    // set the required fields in the ObjectType struct
    PyCustomDustPropertiesType.tp_name = object_name.str().c_str();
    PyCustomDustPropertiesType.tp_basicsize = sizeof(PyCustomDustProperties);
    PyCustomDustPropertiesType.tp_dealloc = (destructor)dealloc;
    PyCustomDustPropertiesType.tp_methods = DustPropertiesMethods;
    PyCustomDustPropertiesType.tp_init = (initproc)init;
    PyCustomDustPropertiesType.tp_new = alloc;

    // finalize creation of the PyCustomDustPropertiesType
    PyType_Ready(&PyCustomDustPropertiesType);
    // add a PyCustomDustPropertiesType to the module
    Py_INCREF(&PyCustomDustPropertiesType);
    PyModule_AddObject(
        module, "CustomDustProperties",
        reinterpret_cast<PyObject *>(&PyCustomDustPropertiesType));
  }
};

#endif // PYDUSTPROPERTIES_HPP
