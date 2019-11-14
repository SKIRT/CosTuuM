/**
 * @file CTMModule.cpp
 *
 * @brief CosTuuM Python module.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "TMatrixCalculator.hpp"

#include <Python.h>
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "structmember.h"

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/// T-matrix object

/**
 * @brief C struct wrapper around a T-matrix object.
 */
typedef struct {
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief T-matrix object. */
  TMatrix *_Tmatrix;
} TmatrixObject;

/**
 * @brief Destructor for the T-matrix object.
 *
 * @param self T-matrix object that is being deallocated.
 */
static void TmatrixObject_dealloc(TmatrixObject *self) {
  // first free the wrapped T-matrix object
  delete self->_Tmatrix;
  // then free the Python object
  Py_TYPE(self)->tp_free((PyObject *)self);
}

/**
 * @brief Constructor for the T-matrix object.
 *
 * Constructs an empty T-matrix object.
 *
 * @param type Type for the object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to the newly created object.
 */
static PyObject *TmatrixObject_new(PyTypeObject *type, PyObject *args,
                                   PyObject *kwargs) {
  // allocate the Python object
  TmatrixObject *self = (TmatrixObject *)type->tp_alloc(type, 0);
  // set the wrapped T-matrix to a sensible initial value
  self->_Tmatrix = nullptr;
  return (PyObject *)self;
}

/**
 * @brief init() function for the T-matrix object.
 *
 * Required arguments are:
 *  - particle_radius: Radius of the scattering particle (in m).
 *  - axis_ratio: Axis ratio of the scattering particle.
 *  - wavelength: Wavelength of the scattered radiation (in m).
 *  - refractive_index: Complex refractive index of the scattering particle.
 *
 * Additional optional arguments are:
 *  - tolerance: Tolerance for the iterative T-matrix calculation.
 *  - maximum_order: Maximum order to use in the spherical basis function
 *    expansion.
 *  - gauss_legendre_factor: Multiplicative factor expressing the number of
 *    Gauss-Legendre quadrature points used during the first phase of the
 *    T-matrix calculation as a multiple of the order for each iteration.
 *  - maximum_ngauss: Maximum number of Gauss-Legendre quadrature points to use
 *    during the T-matrix calculation.
 *
 * @param self T-matrix object that is being initialised.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return 0 on success, 1 on failure.
 */
static int TmatrixObject_init(TmatrixObject *self, PyObject *args,
                              PyObject *kwargs) {

  // required arguments
  float_type particle_radius;
  float_type axis_ratio;
  float_type wavelength;
  std::complex<float_type> mr;

  // optional arguments
  float_type cos2beta = 1. / 3.;
  float_type tolerance = 1.e-4;
  uint_fast32_t maximum_order = 200;
  uint_fast32_t ndgs = 2;
  uint_fast32_t maximum_ngauss = 500;

  // always fixed (for now)
  const float_type ratio_of_radii = 1.;

  /// parse arguments
  // list of keywords (in the expected order)
  // note that we need to use strdup because the Python API expects char*,
  // while C++ strings are const char*
  // not doing this results in compilation warnings
  static char *kwlist[] = {
      strdup("particle_radius"), strdup("axis_ratio"),
      strdup("wavelength"),      strdup("refractive_index"),
      strdup("cos2beta"),        strdup("tolerance"),
      strdup("maximum_order"),   strdup("gauss_legendre_factor"),
      strdup("maximum_ngauss"),  nullptr};

  // allocate a temporary variable to store the complex refractive index
  Py_complex mr_temp;
  // parse the keywords/positional arguments
  // d is a double
  // D is a complex double (Py_complex)
  // I is an unsigned integer
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddD|ddIII", kwlist,
                                   &particle_radius, &axis_ratio, &wavelength,
                                   &mr_temp, &cos2beta, &tolerance,
                                   &maximum_order, &ndgs, &maximum_ngauss)) {
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
  // get the complex components from the Py_complex and store them in the
  // std::complex variable
  mr.real(mr_temp.real);
  mr.imag(mr_temp.imag);

  // construct the T-matrix
  self->_Tmatrix = TMatrixCalculator::calculate_TMatrix(
      ratio_of_radii, axis_ratio, particle_radius, wavelength, maximum_order,
      tolerance, ndgs, mr, maximum_ngauss);

  // average it out over the orientation distribution
  OrientationDistribution *orientation;
  if (cos2beta == 0. || cos2beta == 1.) {
    orientation = new DavisGreensteinOrientationDistribution(
        2 * self->_Tmatrix->get_nmax(), axis_ratio);
  } else if (cos2beta == 1. / 3.) {
    orientation = new OrientationDistribution(2 * self->_Tmatrix->get_nmax());
    orientation->initialise();
  } else {
    orientation = new MishchenkoOrientationDistribution(
        2 * self->_Tmatrix->get_nmax(), axis_ratio, cos2beta);
  }
  self->_Tmatrix = TMatrixCalculator::apply_orientation_distribution(
      *self->_Tmatrix, *orientation);
  delete orientation;

  // everything went well: return 0
  return 0;
}

/**
 * @brief Get the maximum order for a T-matrix object.
 *
 * @param self T-matrix object.
 * @param Py_UNUSED Additional arguments are not used but need to be present.
 * @return Integer Python object containing the maximum order.
 */
static PyObject *TmatrixObject_get_nmax(TmatrixObject *self,
                                        PyObject *Py_UNUSED(ignored)) {
  // wrap the returned value in a Python object
  return PyLong_FromUnsignedLong(self->_Tmatrix->get_nmax());
}

/**
 * @brief Compute the extinction matrix using the given T-matrix object.
 *
 * Required arguments are:
 *  - theta: Input zenith angle (in radians).
 *  - phi: Input azimuth angle (in radians).
 *
 * Additional optional arguments are:
 *  - alpha: Azimuth angle of the particle rotation axis (in radians).
 *  - beta: Zenith angle of the particle rotation axis (in radians).
 *
 * @param self T-matrix object being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to a 4x4 double precision NumPy array containing the
 * components of the extinction matrix.
 */
static PyObject *TmatrixObject_get_extinction_matrix(TmatrixObject *self,
                                                     PyObject *args,
                                                     PyObject *kwargs) {

  // required arguments
  float_type theta;
  float_type phi;

  // optional arguments
  float_type alpha = 0.;
  float_type beta = 0.;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("theta"), strdup("phi"), strdup("alpha"),
                           strdup("beta"), nullptr};

  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|dd", kwlist, &theta, &phi,
                                   &alpha, &beta)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }

  // get the extinction matrix
  Matrix<float_type> K =
      self->_Tmatrix->get_extinction_matrix(alpha, beta, theta, phi);

  // create an uninitialised 4x4 double NumPy array
  npy_intp dims[2] = {4, 4};
  PyArrayObject *numpy_array =
      (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  // copy the elements from the extinction matrix into the array
  for (uint_fast8_t i = 0; i < 4; ++i) {
    for (uint_fast8_t j = 0; j < 4; ++j) {
      *((float_type *)PyArray_GETPTR2(numpy_array, i, j)) = K(i, j);
    }
  }

  // return the array
  return PyArray_Return(numpy_array);
}

/*! @brief List of T-matrix object methods that are exposed to Python. */
static PyMethodDef TmatrixObject_methods[] = {
    {"get_nmax", (PyCFunction)TmatrixObject_get_nmax, METH_NOARGS,
     "Return the maximum order of the T-matrix."},
    {"get_extinction_matrix", (PyCFunction)TmatrixObject_get_extinction_matrix,
     METH_VARARGS | METH_KEYWORDS, "Return the extinction matrix."},
    {nullptr}};

/*! @brief Python Object type for the T-matrix (is edited in the module
 *  initialisation function since C++ does not allow fancy unordered struct
 *  initialisation). */
static PyTypeObject TmatrixObjectType = {PyVarObject_HEAD_INIT(nullptr, 0)};

/// MishchenkoOrientationDistributionObject

/**
 * @brief C struct wrapper around a MishchenkoOrientationDistribution object.
 */
typedef struct {
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief T-matrix object. */
  MishchenkoOrientationDistribution *_distribution;
} MishchenkoOrientationDistributionObject;

/**
 * @brief Destructor for the MishchenkoOrientationDistributionObject.
 *
 * @param self MishchenkoOrientationDistributionObject that is being
 * deallocated.
 */
static void MishchenkoOrientationDistributionObject_dealloc(
    MishchenkoOrientationDistributionObject *self) {
  // first free the wrapped T-matrix object
  delete self->_distribution;
  // then free the Python object
  Py_TYPE(self)->tp_free((PyObject *)self);
}

/**
 * @brief Constructor for the MishchenkoOrientationDistributionObject.
 *
 * Constructs an empty MishchenkoOrientationDistributionObject.
 *
 * @param type Type for the object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to the newly created object.
 */
static PyObject *MishchenkoOrientationDistributionObject_new(PyTypeObject *type,
                                                             PyObject *args,
                                                             PyObject *kwargs) {
  // allocate the Python object
  MishchenkoOrientationDistributionObject *self =
      (MishchenkoOrientationDistributionObject *)type->tp_alloc(type, 0);
  // set the wrapped T-matrix to a sensible initial value
  self->_distribution = nullptr;
  return (PyObject *)self;
}

/**
 * @brief init() function for the MishchenkoOrientationDistributionObject.
 *
 * Required arguments are:
 *  - cos2beta: Parameter for the object.
 *
 * Additional optional arguments are:
 *  - maximum_order: Maximum order to use in the spherical basis function
 *    expansion.
 *
 * @param self MishchenkoOrientationDistributionObject that is being
 * initialised.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return 0 on success, 1 on failure.
 */
static int MishchenkoOrientationDistributionObject_init(
    MishchenkoOrientationDistributionObject *self, PyObject *args,
    PyObject *kwargs) {

  // required arguments
  float_type axis_ratio;
  float_type cos2beta;

  // optional arguments
  uint_fast32_t maximum_order = 200;

  /// parse arguments
  // list of keywords (in the expected order)
  // note that we need to use strdup because the Python API expects char*,
  // while C++ strings are const char*
  // not doing this results in compilation warnings
  static char *kwlist[] = {strdup("axis_ratio"), strdup("cos2beta"),
                           strdup("maximum_order"), nullptr};

  // parse the keywords/positional arguments
  // d is a double
  // I is an unsigned integer
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|I", kwlist, &axis_ratio,
                                   &cos2beta, &maximum_order)) {
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

  self->_distribution = new MishchenkoOrientationDistribution(
      maximum_order, axis_ratio, cos2beta);

  // everything went well: return 0
  return 0;
}

/**
 * @brief Compute the orientation distribution for the given angle.
 *
 * Required arguments are:
 *  - beta: Zenith angle of the particle rotation axis (in radians).
 *
 * @param self MishchenkoOrientationDistributionObject being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Value of the distribution, wrapped in a Python float object.
 */
static PyObject *MishchenkoOrientationDistributionObject_get_distribution(
    MishchenkoOrientationDistributionObject *self, PyObject *args,
    PyObject *kwargs) {

  // required arguments
  float_type beta;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("beta"), nullptr};

  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &beta)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }

  const float_type value =
      (*static_cast<OrientationDistribution *>(self->_distribution))(beta);

  // return the array
  return PyFloat_FromDouble(value);
}

/*! @brief List of MishchenkoOrientationDistributionObject methods that are
 *  exposed to Python. */
static PyMethodDef MishchenkoOrientationDistributionObject_methods[] = {
    {"get_distribution",
     (PyCFunction)MishchenkoOrientationDistributionObject_get_distribution,
     METH_VARARGS | METH_KEYWORDS, "Return the distribution function."},
    {nullptr}};

/*! @brief Python Object type for the MishchenkoOrientationDistributionObject
 *  (is edited in the module initialisation function since C++ does not allow
 *  fancy unordered struct initialisation). */
static PyTypeObject MishchenkoOrientationDistributionObjectType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/// CosTuuM module

/*! @brief Module definition for the CTM module. */
static struct PyModuleDef CTMmodule = {PyModuleDef_HEAD_INIT, "CosTuuM",
                                       "T-matrix module.", -1, 0};

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

  // set the required fields in the TmatrixObjectType struct
  // we need to do this because C++ does not support fancy C struct
  // initialisation
  // without this, we would need to manually set all elements in the struct
  // this is a pain (there are many) and makes it very hard to make the code
  // portable (addition of a single element to the struct would break it)
  TmatrixObjectType.tp_name = "CosTuuM.TMatrix";
  TmatrixObjectType.tp_basicsize = sizeof(TmatrixObject);
  TmatrixObjectType.tp_dealloc = (destructor)TmatrixObject_dealloc;
  TmatrixObjectType.tp_methods = TmatrixObject_methods;
  TmatrixObjectType.tp_init = (initproc)TmatrixObject_init;
  TmatrixObjectType.tp_new = TmatrixObject_new;

  // finalize creation of the TmatrixObjectType
  PyType_Ready(&TmatrixObjectType);
  // add a T-matrix object to the module
  // any call to CTMmodule.TMatrix() will use this same object
  Py_INCREF(&TmatrixObjectType);
  PyModule_AddObject(m, "TMatrix", (PyObject *)&TmatrixObjectType);

  // set the required fields in the MishchenkoOrientationDistributionObjectType
  // struct
  MishchenkoOrientationDistributionObjectType.tp_name =
      "CosTuuM.MishchenkoOrientationDistribution";
  MishchenkoOrientationDistributionObjectType.tp_basicsize =
      sizeof(MishchenkoOrientationDistributionObject);
  MishchenkoOrientationDistributionObjectType.tp_dealloc =
      (destructor)MishchenkoOrientationDistributionObject_dealloc;
  MishchenkoOrientationDistributionObjectType.tp_methods =
      MishchenkoOrientationDistributionObject_methods;
  MishchenkoOrientationDistributionObjectType.tp_init =
      (initproc)MishchenkoOrientationDistributionObject_init;
  MishchenkoOrientationDistributionObjectType.tp_new =
      MishchenkoOrientationDistributionObject_new;

  // finalize creation of the MishchenkoOrientationDistributionObjectType
  PyType_Ready(&MishchenkoOrientationDistributionObjectType);
  // add a MishchenkoOrientationDistributionObject to the module
  // any call to CTMmodule.MishchenkoOrientationDistribution() will use this
  // same object
  Py_INCREF(&MishchenkoOrientationDistributionObjectType);
  PyModule_AddObject(m, "MishchenkoOrientationDistribution",
                     (PyObject *)&MishchenkoOrientationDistributionObjectType);

  // return the module object
  return m;
}
