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
#include <sstream>

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
   * @param self DustProperties wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyDustProperties *self, PyObject *args, PyObject *kwargs) {

    self->_dust_properties = new DraineDustProperties();

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

public:
  /**
   * @brief Constructor.
   *
   * @param python_function Python function that is being wrapped.
   */
  inline CustomDustProperties(PyObject *python_function)
      : _python_function(python_function) {
    // tell Python we store a copy of the function object
    Py_INCREF(_python_function);
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~CustomDustProperties() {
    // tell Python we are done with our local copy of the function object
    Py_DECREF(_python_function);
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
    double wavelength_d = wavelength;
    double grain_size_d = grain_size;
    // wrappers for int values
    int grain_type_i = grain_type;
    // build the argument list for the Python function call
    PyObject *arglist =
        Py_BuildValue("(ddi)", &wavelength_d, &grain_size_d, &grain_type_i);
    // place the Python function call
    PyObject *result = PyEval_CallObject(_python_function, arglist);
    // we are done with the argument list
    Py_DECREF(arglist);

    // parse the result
    const double result_real = PyComplex_RealAsDouble(result);
    const double result_imag = PyComplex_ImagAsDouble(result);
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
   * @param self DustProperties wrapper object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyDustProperties *self, PyObject *args, PyObject *kwargs) {

    // required arguments
    PyObject *function;

    // list of keywords
    static char *kwlist[] = {strdup("function"), nullptr};

    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist, &function)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      return 1;
    }

    self->_dust_properties = new CustomDustProperties(function);

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
