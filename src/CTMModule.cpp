/**
 * @file CTMModule.cpp
 *
 * @brief CosTuuM Python module.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
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
  delete self->_Tmatrix;
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
  TmatrixObject *self = (TmatrixObject *)type->tp_alloc(type, 0);
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

  float_type particle_radius;
  float_type axis_ratio;
  float_type wavelength;
  std::complex<float_type> mr;

  float_type tolerance = 1.e-4;
  uint_fast32_t maximum_order = 200;
  uint_fast32_t ndgs = 2;
  uint_fast32_t maximum_ngauss = 500;

  const float_type ratio_of_radii = 1.;

  static char *kwlist[] = {strdup("particle_radius"),
                           strdup("axis_ratio"),
                           strdup("wavelength"),
                           strdup("refractive_index"),
                           strdup("tolerance"),
                           strdup("maximum_order"),
                           strdup("gauss_legendre_factor"),
                           strdup("maximum_ngauss"),
                           nullptr};

  Py_complex mr_temp;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddD|dIII", kwlist,
                                   &particle_radius, &axis_ratio, &wavelength,
                                   &mr_temp, &tolerance, &maximum_order, &ndgs,
                                   &maximum_ngauss)) {
    ctm_warning("Wrong arguments provided!");
    return 1;
  }
  mr.real(mr_temp.real);
  mr.imag(mr_temp.imag);

  self->_Tmatrix = TMatrixCalculator::calculate_TMatrix(
      ratio_of_radii, axis_ratio, particle_radius, wavelength, maximum_order,
      tolerance, ndgs, mr, maximum_ngauss);

  OrientationDistribution orientation(2 * self->_Tmatrix->get_nmax());
  self->_Tmatrix = TMatrixCalculator::apply_orientation_distribution(
      *self->_Tmatrix, orientation);

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

  float_type theta;
  float_type phi;

  float_type alpha = 0.;
  float_type beta = 0.;

  static char *kwlist[] = {strdup("theta"), strdup("phi"), strdup("alpha"),
                           strdup("beta"), nullptr};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd|dd", kwlist, &theta, &phi,
                                   &alpha, &beta)) {
    ctm_warning("Wrong arguments provided!");
    return nullptr;
  }

  Matrix<float_type> K =
      self->_Tmatrix->get_extinction_matrix(alpha, beta, theta, phi);

  npy_intp dims[2] = {4, 4};
  PyArrayObject *numpy_array =
      (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  for (uint_fast8_t i = 0; i < 4; ++i) {
    for (uint_fast8_t j = 0; j < 4; ++j) {
      *((float_type *)PyArray_GETPTR2(numpy_array, i, j)) = K(i, j);
    }
  }

  return PyArray_Return(numpy_array);
}

/*! @brief List of T-matrix object methods that are exposed to Python. */
static PyMethodDef TmatrixObject_methods[] = {
    {"get_nmax", (PyCFunction)TmatrixObject_get_nmax, METH_NOARGS,
     "Return the maximum order of the T-matrix."},
    {"get_extinction_matrix", (PyCFunction)TmatrixObject_get_extinction_matrix,
     METH_VARARGS | METH_KEYWORDS, "Return the excintion matrix."},
    {nullptr}};

/*! @brief Python Object type for the T-matrix (is edited in the module
 *  initialisation function since C++ does not allow fancy unordered struct
 *  initialisation). */
static PyTypeObject TmatrixObjectType = {PyVarObject_HEAD_INIT(nullptr, 0)};

/*! @brief Module definition for the CTM module. */
static struct PyModuleDef CTMmodule = {PyModuleDef_HEAD_INIT, "CTMmodule",
                                       "T-matrix module.", -1, 0};

/**
 * @brief CTMmodule initialisation function.
 *
 * @return Pointer to the initialised module object.
 */
PyMODINIT_FUNC PyInit_CTMmodule() {

  import_array();
  PyObject *m = PyModule_Create(&CTMmodule);

  TmatrixObjectType.tp_name = "CTMmodule.TMatrix";
  TmatrixObjectType.tp_basicsize = sizeof(TmatrixObject);
  TmatrixObjectType.tp_dealloc = (destructor)TmatrixObject_dealloc;
  TmatrixObjectType.tp_methods = TmatrixObject_methods;
  TmatrixObjectType.tp_init = (initproc)TmatrixObject_init;
  TmatrixObjectType.tp_new = TmatrixObject_new;

  PyType_Ready(&TmatrixObjectType);
  Py_INCREF(&TmatrixObjectType);
  PyModule_AddObject(m, "TMatrix", (PyObject *)&TmatrixObjectType);
  return m;
}
