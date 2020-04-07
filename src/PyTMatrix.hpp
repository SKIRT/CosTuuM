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
 * @file PyTMatrix.hpp
 *
 * @brief Python exposure of the (old) TMatrix object.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYTMATRIX_HPP
#define PYTMATRIX_HPP

#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "PyHelperFunctions.hpp"
#include "TMatrixCalculator.hpp"

#include <Python.h>
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/*! @brief Python Object type for the T-matrix (is edited in the module
 *  initialisation function since C++ does not allow fancy unordered struct
 *  initialisation). */
static PyTypeObject TmatrixObjectType = {PyVarObject_HEAD_INIT(nullptr, 0)};

/**
 * @brief Python wrapper around a T-matrix object.
 */
class PyTMatrix {
public:
  /*! @brief Python object members. */
  PyObject_HEAD;

  /*! @brief T-matrix object. */
  TMatrix *_Tmatrix;

  /**
   * @brief Destructor for the T-matrix object.
   *
   * @param self T-matrix object that is being deallocated.
   */
  static void dealloc(PyTMatrix *self) {
    // first free the wrapped T-matrix object
    delete self->_Tmatrix;
    // then free the Python object
    Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
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
  static PyObject *alloc(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    // allocate the Python object
    PyTMatrix *self = reinterpret_cast<PyTMatrix *>(type->tp_alloc(type, 0));
    // set the wrapped T-matrix to a sensible initial value
    self->_Tmatrix = nullptr;
    return reinterpret_cast<PyObject *>(self);
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
   *  - maximum_ngauss: Maximum number of Gauss-Legendre quadrature points to
   * use during the T-matrix calculation.
   *  - is_equal_volume_radius: Is the input particle radius an equal volume
   *    radius?
   *
   * @param self T-matrix object that is being initialised.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return 0 on success, 1 on failure.
   */
  static int init(PyTMatrix *self, PyObject *args, PyObject *kwargs) {

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
    float_type ratio_of_radii = 1.;

    /// parse arguments
    // list of keywords (in the expected order)
    // note that we need to use strdup because the Python API expects char*,
    // while C++ strings are const char*
    // not doing this results in compilation warnings
    static char *kwlist[] = {strdup("particle_radius"),
                             strdup("axis_ratio"),
                             strdup("wavelength"),
                             strdup("refractive_index"),
                             strdup("cos2beta"),
                             strdup("tolerance"),
                             strdup("maximum_order"),
                             strdup("gauss_legendre_factor"),
                             strdup("maximum_ngauss"),
                             strdup("is_equal_volume_radius"),
                             nullptr};

    // allocate temporary variables to store double precision arguments
    double particle_radius_d, axis_ratio_d, wavelength_d;
    double cos2beta_d = static_cast<double>(cos2beta);
    double tolerance_d = static_cast<double>(tolerance);
    // temporary variables to store integer arguments
    unsigned int maximum_order_i = maximum_order;
    unsigned int ndgs_i = ndgs;
    unsigned int maximum_ngauss_i = maximum_ngauss;
    // allocate a temporary variable to store the complex refractive index
    Py_complex mr_temp;
    // allocate a temporary variable to store the equal volume radius boolean
    int is_equal_volume_radius = true;
    // parse the keywords/positional arguments
    // d is a double
    // D is a complex double (Py_complex)
    // I is an unsigned integer
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "dddD|ddIIIp", kwlist, &particle_radius_d,
            &axis_ratio_d, &wavelength_d, &mr_temp, &cos2beta_d, &tolerance_d,
            &maximum_order_i, &ndgs_i, &maximum_ngauss_i,
            &is_equal_volume_radius)) {
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
    particle_radius = particle_radius_d;
    axis_ratio = axis_ratio_d;
    wavelength = wavelength_d;
    cos2beta = cos2beta_d;
    tolerance = tolerance_d;
    // unpack integer arguments
    maximum_order = maximum_order_i;
    ndgs = ndgs_i;
    maximum_ngauss = maximum_ngauss_i;
    // get the complex components from the Py_complex and store them in the
    // std::complex variable
    mr.real(mr_temp.real);
    mr.imag(mr_temp.imag);
    // set the ratio_of_radii if it is not correct
    if (!is_equal_volume_radius) {
      ratio_of_radii = 0.1;
    }

    // construct the T-matrix
    self->_Tmatrix = TMatrixCalculator::calculate_TMatrix(
        ratio_of_radii, axis_ratio, particle_radius, wavelength, maximum_order,
        tolerance, ndgs, mr, maximum_ngauss);

    // average it out over the orientation distribution
    if (cos2beta >= 0.) {
      OrientationDistribution *orientation;
      if (cos2beta == 0. || cos2beta == 1.) {
        orientation = new DavisGreensteinOrientationDistribution(
            2 * self->_Tmatrix->get_nmax(), axis_ratio);
      } else if (cos2beta == 1. / 3.) {
        orientation =
            new OrientationDistribution(2 * self->_Tmatrix->get_nmax());
        orientation->initialise();
      } else {
        orientation = new MishchenkoOrientationDistribution(
            2 * self->_Tmatrix->get_nmax(), cos2beta);
      }
      self->_Tmatrix = TMatrixCalculator::apply_orientation_distribution(
          *self->_Tmatrix, *orientation);
      delete orientation;
    }

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
  static PyObject *get_nmax(PyTMatrix *self, PyObject *Py_UNUSED(ignored)) {
    // wrap the returned value in a Python object
    return PyLong_FromUnsignedLong(self->_Tmatrix->get_nmax());
  }

  /**
   * @brief Get the absorption cross section appropriately averaged over the
   * outgoing angle @f$\theta{}@f$.
   *
   * @param self T-matrix object.
   * @param Py_UNUSED Additional arguments are not used but need to be present.
   * @return Integer Python object containing the maximum order.
   */
  static PyObject *
  get_average_absorption_cross_section(PyTMatrix *self,
                                       PyObject *Py_UNUSED(ignored)) {
    // wrap the returned value in a Python object
    return PyFloat_FromDouble(static_cast<double>(
        self->_Tmatrix->get_average_absorption_cross_section(20)));
  }

  /**
   * @brief Get the absorption cross section for the given input angle(s).
   *
   * Required arguments are:
   *  - theta: Input zenith angle (in radians).
   *
   * Additional optional arguments are:
   *  - ngauss: Number of Gauss-Legendre quadrature points to use to compute the
   *    integral to subtract the scattering contribution to the extinction.
   *
   * @param self T-matrix object being used.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Pointer to an absorption cross section value or array that has the
   * same shape as the input angles.
   */
  static PyObject *get_absorption_cross_section(PyTMatrix *self, PyObject *args,
                                                PyObject *kwargs) {

    // required arguments
    PyArrayObject *thetas;

    // optional arguments
    uint_fast32_t ngauss = 100;

    // list of keywords (see comment above)
    static char *kwlist[] = {strdup("theta"), strdup("ngauss"), nullptr};

    // temporary variables to store integer arguments
    unsigned int ngauss_i = ngauss;
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|I", kwlist,
                                     PyArray_Converter, &thetas, &ngauss_i)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }
    // unpack integer variables
    ngauss = ngauss_i;

    // determine the size of the input arrays
    npy_intp thetasize;

    // get the number of dimensions for the theta array
    const npy_intp thetadim = PyArray_NDIM(thetas);
    // we only accept 0D (scalar) or 1D input
    if (thetadim > 1) {
      ctm_warning("Wrong shape for input array!");
      return nullptr;
    }
    // check if we are in the 0D or 1D case
    if (thetadim > 0) {
      // if 1D, simply get the size from the 1 dimension
      const npy_intp *thetadims = PyArray_DIMS(thetas);
      thetasize = thetadims[0];
    } else {
      // if 0D, set the size to 1
      thetasize = 1;
      // reshape the array into a 1D array with 1 element, so that we can
      // manipulate it in the same way as a 1D array
      npy_intp newdims[1] = {1};
      PyArray_Dims newdimsobj;
      newdimsobj.ptr = newdims;
      newdimsobj.len = 1;
      thetas = reinterpret_cast<PyArrayObject *>(
          PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
    }

    // create an uninitialised double NumPy array to store the results
    // shape: thetasize
    npy_intp dims[1] = {thetasize};
    PyArrayObject *Cabs = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(1, dims, NPY_DOUBLE));

    const double *thetadata = reinterpret_cast<double *>(PyArray_DATA(thetas));
    std::vector<float_type> thetadata_f(thetasize);
    for (npy_intp i = 0; i < thetasize; ++i) {
      thetadata_f[i] = thetadata[i];
    }

    double *Cabsdata = reinterpret_cast<double *>(PyArray_DATA(Cabs));
    std::vector<float_type> Cabsdata_f(thetasize);
    self->_Tmatrix->get_absorption_cross_section(&thetadata_f[0], thetasize,
                                                 ngauss, &Cabsdata_f[0]);
    for (npy_intp i = 0; i < thetasize; ++i) {
      Cabsdata[i] = static_cast<double>(Cabsdata_f[i]);
    }

    // tell Python we are done with the input objects (so that memory is
    // properly deallocated)
    Py_DECREF(thetas);

    return PyArray_Return(Cabs);
  }

  /**
   * @brief Get the absorption cross sections for the given input angle(s).
   *
   * Required arguments are:
   *  - theta: Input zenith angle (in radians).
   *
   * Additional optional arguments are:
   *  - ngauss: Number of Gauss-Legendre quadrature points to use to compute the
   *    integral to subtract the scattering contribution to the extinction.
   *
   * @param self T-matrix object being used.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Pointer to an array containing the cross section values. The array
   * has two elements per input angle theta.
   */
  static PyObject *get_absorption_cross_sections(PyTMatrix *self,
                                                 PyObject *args,
                                                 PyObject *kwargs) {

    // required arguments
    PyArrayObject *thetas;

    // optional arguments
    uint_fast32_t ngauss = 100;

    // list of keywords (see comment above)
    static char *kwlist[] = {strdup("theta"), strdup("ngauss"), nullptr};

    // temporary variables to store integer arguments
    unsigned int ngauss_i = ngauss;
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|I", kwlist,
                                     PyArray_Converter, &thetas, &ngauss_i)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }
    // unpack integer arguments
    ngauss = ngauss_i;

    // determine the size of the input arrays
    npy_intp thetasize;

    // get the number of dimensions for the theta array
    const npy_intp thetadim = PyArray_NDIM(thetas);
    // we only accept 0D (scalar) or 1D input
    if (thetadim > 1) {
      ctm_warning("Wrong shape for input array!");
      return nullptr;
    }
    // check if we are in the 0D or 1D case
    if (thetadim > 0) {
      // if 1D, simply get the size from the 1 dimension
      const npy_intp *thetadims = PyArray_DIMS(thetas);
      thetasize = thetadims[0];
    } else {
      // if 0D, set the size to 1
      thetasize = 1;
      // reshape the array into a 1D array with 1 element, so that we can
      // manipulate it in the same way as a 1D array
      npy_intp newdims[1] = {1};
      PyArray_Dims newdimsobj;
      newdimsobj.ptr = newdims;
      newdimsobj.len = 1;
      thetas = reinterpret_cast<PyArrayObject *>(
          PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
    }

    // create an uninitialised double NumPy array to store the results
    // shape: thetasize x 2
    npy_intp dims[2] = {thetasize, 2};
    PyArrayObject *Cabs = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(2, dims, NPY_DOUBLE));

    const double *thetadata = reinterpret_cast<double *>(PyArray_DATA(thetas));
    std::vector<float_type> thetadata_f(thetasize);
    for (npy_intp i = 0; i < thetasize; ++i) {
      thetadata_f[i] = thetadata[i];
    }
    double *Cabsdata = reinterpret_cast<double *>(PyArray_DATA(Cabs));
    std::vector<float_type> Cabsdata_f(2 * thetasize);
    self->_Tmatrix->get_absorption_cross_sections(&thetadata_f[0], thetasize,
                                                  ngauss, &Cabsdata_f[0]);
    for (npy_intp i = 0; i < 2 * thetasize; ++i) {
      Cabsdata[i] = static_cast<double>(Cabsdata_f[i]);
    }

    // tell Python we are done with the input objects (so that memory is
    // properly deallocated)
    Py_DECREF(thetas);

    return PyArray_Squeeze(Cabs);
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
  static PyObject *get_extinction_matrix(PyTMatrix *self, PyObject *args,
                                         PyObject *kwargs) {

    // required arguments
    PyArrayObject *thetas;
    PyArrayObject *phis;

    // optional arguments
    float_type alpha = 0.;
    float_type beta = 0.;

    // list of keywords (see comment above)
    static char *kwlist[] = {strdup("theta"), strdup("phi"), strdup("alpha"),
                             strdup("beta"), nullptr};

    // allocate temporary variables to store double precision arguments
    double alpha_d = static_cast<double>(alpha);
    double beta_d = static_cast<double>(beta);
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O&O&|dd", kwlist, PyArray_Converter, &thetas,
            PyArray_Converter, &phis, &alpha_d, &beta_d)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }
    // unpack double precision values
    alpha = alpha_d;
    beta = beta_d;

    // determine the size of the input arrays
    npy_intp thetasize;
    npy_intp phisize;

    // get the number of dimensions for the theta array
    const npy_intp thetadim = PyArray_NDIM(thetas);
    // we only accept 0D (scalar) or 1D input
    if (thetadim > 1) {
      ctm_warning("Wrong shape for input array!");
      return nullptr;
    }
    // check if we are in the 0D or 1D case
    if (thetadim > 0) {
      // if 1D, simply get the size from the 1 dimension
      const npy_intp *thetadims = PyArray_DIMS(thetas);
      thetasize = thetadims[0];
    } else {
      // if 0D, set the size to 1
      thetasize = 1;
      // reshape the array into a 1D array with 1 element, so that we can
      // manipulate it in the same way as a 1D array
      npy_intp newdims[1] = {1};
      PyArray_Dims newdimsobj;
      newdimsobj.ptr = newdims;
      newdimsobj.len = 1;
      thetas = reinterpret_cast<PyArrayObject *>(
          PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
    }
    // now repeat for the phi array
    const npy_intp phidim = PyArray_NDIM(phis);
    if (phidim > 1) {
      ctm_warning("Wrong shape for input array!");
      return nullptr;
    }
    if (phidim > 0) {
      const npy_intp *phidims = PyArray_DIMS(phis);
      phisize = phidims[0];
    } else {
      phisize = 1;
      npy_intp newdims[1] = {1};
      PyArray_Dims newdimsobj;
      newdimsobj.ptr = newdims;
      newdimsobj.len = 1;
      phis = reinterpret_cast<PyArrayObject *>(
          PyArray_Newshape(phis, &newdimsobj, NPY_ANYORDER));
    }

    // create an uninitialised double NumPy array to store the results
    // shape: thetasize x phisize x 4 x 4
    npy_intp dims[4] = {thetasize, phisize, 4, 4};
    PyArrayObject *numpy_array = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(4, dims, NPY_DOUBLE));

    // loop over all theta elements
    for (npy_intp itheta = 0; itheta < thetasize; ++itheta) {
      // get the corresponding theta angle
      const float_type theta =
          *(reinterpret_cast<double *>(PyArray_GETPTR1(thetas, itheta)));
      // loop over all phi elements
      for (npy_intp iphi = 0; iphi < phisize; ++iphi) {
        // get the corresponding phi angle
        const float_type phi =
            *(reinterpret_cast<double *>(PyArray_GETPTR1(phis, iphi)));

        // get the extinction matrix
        Matrix<float_type> K =
            self->_Tmatrix->get_extinction_matrix(alpha, beta, theta, phi);

        // copy the elements from the extinction matrix into the array
        for (uint_fast8_t i = 0; i < 4; ++i) {
          for (uint_fast8_t j = 0; j < 4; ++j) {
            *reinterpret_cast<double *>(
                PyArray_GETPTR4(numpy_array, itheta, iphi, i, j)) =
                static_cast<double>(K(i, j));
          }
        }
      }
    }

    // tell Python we are done with the input objects (so that memory is
    // properly deallocated)
    Py_DECREF(thetas);
    Py_DECREF(phis);

    // return the array
    // we squeeze it so that the output dimensions match the input arrays:
    // a 1D theta and 0D phi array e.g. will have a thetasize x 4 x 4 output
    return PyArray_Squeeze(numpy_array);
  }

  /**
   * @brief Compute the scattering matrix using the given T-matrix object.
   *
   * Required arguments are:
   *  - theta_in: Input zenith angle (in radians).
   *  - phi_in: Input azimuth angle (in radians).
   *  - theta_out: Output zenith angle (in radians).
   *  - phi_out: Output azimuth angle (in radians).
   *
   * Additional optional arguments are:
   *  - alpha: Azimuth angle of the particle rotation axis (in radians).
   *  - beta: Zenith angle of the particle rotation axis (in radians).
   *
   * @param self T-matrix object being used.
   * @param args Positional arguments.
   * @param kwargs Keyword arguments.
   * @return Pointer to a 4x4 double precision NumPy array containing the
   * components of the scattering matrix.
   */
  static PyObject *get_scattering_matrix(PyTMatrix *self, PyObject *args,
                                         PyObject *kwargs) {

    // required arguments
    PyArrayObject *pytheta_ins;
    PyArrayObject *pyphi_ins;
    PyArrayObject *pytheta_outs;
    PyArrayObject *pyphi_outs;

    // optional arguments
    float_type alpha = 0.;
    float_type beta = 0.;

    // list of keywords (see comment above)
    static char *kwlist[] = {strdup("theta_in"),
                             strdup("phi_in"),
                             strdup("theta_out"),
                             strdup("phi_out"),
                             strdup("alpha"),
                             strdup("beta"),
                             nullptr};

    // allocate temporary variables to store double precision arguments
    double alpha_d = static_cast<double>(alpha);
    double beta_d = static_cast<double>(beta);
    // parse positional and keyword arguments
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O&O&O&O&|dd", kwlist, PyArray_Converter,
            &pytheta_ins, PyArray_Converter, &pyphi_ins, PyArray_Converter,
            &pytheta_outs, PyArray_Converter, &pyphi_outs, &alpha_d, &beta_d)) {
      // again, we do not call ctm_error to avoid killing the Python interpreter
      ctm_warning("Wrong arguments provided!");
      // this time, a nullptr return will signal an error to Python
      return nullptr;
    }
    // unpack double precision values
    alpha = alpha_d;
    beta = beta_d;

    // convert the input arrays to C++ vectors
    std::vector<float_type> theta_ins =
        unpack_numpy_array<float_type, double>(pytheta_ins);
    Py_DECREF(pytheta_ins);
    std::vector<float_type> phi_ins =
        unpack_numpy_array<float_type, double>(pyphi_ins);
    Py_DECREF(pyphi_ins);
    std::vector<float_type> theta_outs =
        unpack_numpy_array<float_type, double>(pytheta_outs);
    Py_DECREF(pytheta_outs);
    std::vector<float_type> phi_outs =
        unpack_numpy_array<float_type, double>(pyphi_outs);
    Py_DECREF(pyphi_outs);

    const npy_intp theta_in_size = theta_ins.size();
    const npy_intp phi_in_size = phi_ins.size();
    const npy_intp theta_out_size = theta_outs.size();
    const npy_intp phi_out_size = phi_outs.size();

    // create an uninitialised double NumPy array to store the results
    // shape: theta_in_size x phi_in_size x theta_out_size x phi_out_size x 4 x
    // 4
    npy_intp dims[6] = {theta_in_size, phi_in_size, theta_out_size,
                        phi_out_size,  4,           4};
    PyArrayObject *numpy_array = reinterpret_cast<PyArrayObject *>(
        PyArray_SimpleNew(6, dims, NPY_DOUBLE));

    // loop over all theta_in elements
    for (npy_intp itheta_in = 0; itheta_in < theta_in_size; ++itheta_in) {
      // get the corresponding theta_in angle
      const float_type theta_in = theta_ins[itheta_in];
      // loop over all phi_in elements
      for (npy_intp iphi_in = 0; iphi_in < phi_in_size; ++iphi_in) {
        // get the corresponding phi_in angle
        const float_type phi_in = phi_ins[iphi_in];

        // loop over all theta_out elements
        for (npy_intp itheta_out = 0; itheta_out < theta_out_size;
             ++itheta_out) {
          // get the corresponding theta_out angle
          const float_type theta_out = theta_outs[itheta_out];
          // loop over all phi_out elements
          for (npy_intp iphi_out = 0; iphi_out < phi_out_size; ++iphi_out) {
            // get the corresponding phi_out angle
            const float_type phi_out = phi_outs[iphi_out];

            // get the scattering matrix
            Matrix<float_type> Z = self->_Tmatrix->get_scattering_matrix(
                alpha, beta, theta_in, phi_in, theta_out, phi_out);

            // copy the elements from the extinction matrix into the array
            for (uint_fast8_t irow = 0; irow < 4; ++irow) {
              for (uint_fast8_t icol = 0; icol < 4; ++icol) {
                npy_intp array_index[6] = {itheta_in, iphi_in, itheta_out,
                                           iphi_out,  irow,    icol};
                *reinterpret_cast<double *>(
                    PyArray_GetPtr(numpy_array, array_index)) =
                    static_cast<double>(Z(irow, icol));
              }
            }
          }
        }
      }
    }

    // return the array
    // we squeeze it so that the output dimensions match the input arrays:
    // a 1D theta and 0D phi array e.g. will have a thetasize x 4 x 4 output
    return PyArray_Squeeze(numpy_array);
  }

  /**
   * @brief Initialize the T-matrix object and add it to the given module.
   *
   * @param module Module to add the object to.
   */
  inline static void initialize(PyObject *module);
};

/*! @brief List of T-matrix object methods that are exposed to Python. */
static PyMethodDef TmatrixObject_methods[] = {
    {"get_nmax", reinterpret_cast<PyCFunction>(PyTMatrix::get_nmax),
     METH_NOARGS, "Return the maximum order of the T-matrix."},
    {"get_extinction_matrix",
     reinterpret_cast<PyCFunction>(PyTMatrix::get_extinction_matrix),
     METH_VARARGS | METH_KEYWORDS, "Return the extinction matrix."},
    {"get_average_absorption_cross_section",
     reinterpret_cast<PyCFunction>(
         PyTMatrix::get_average_absorption_cross_section),
     METH_NOARGS,
     "Return the angular average of the absorption cross section."},
    {"get_absorption_cross_section",
     reinterpret_cast<PyCFunction>(PyTMatrix::get_absorption_cross_section),
     METH_VARARGS | METH_KEYWORDS,
     "Return the absorption cross section for the given input angle(s)."},
    {"get_absorption_cross_sections",
     reinterpret_cast<PyCFunction>(PyTMatrix::get_absorption_cross_sections),
     METH_VARARGS | METH_KEYWORDS,
     "Return both absorption cross sections for the given input angle(s)."},
    {"get_scattering_matrix",
     reinterpret_cast<PyCFunction>(PyTMatrix::get_scattering_matrix),
     METH_VARARGS | METH_KEYWORDS,
     "Return the scattering matrix for the given scattering event."},
    {nullptr}};

void PyTMatrix::initialize(PyObject *module) {

  std::stringstream object_name;
  object_name << PyModule_GetName(module);
  object_name << ".TMatrix";
  // set the required fields in the TmatrixObjectType struct
  // we need to do this because C++ does not support fancy C struct
  // initialisation
  // without this, we would need to manually set all elements in the struct
  // this is a pain (there are many) and makes it very hard to make the code
  // portable (addition of a single element to the struct would break it)
  TmatrixObjectType.tp_name = object_name.str().c_str();
  TmatrixObjectType.tp_basicsize = sizeof(PyTMatrix);
  TmatrixObjectType.tp_dealloc = (destructor)dealloc;
  TmatrixObjectType.tp_methods = TmatrixObject_methods;
  TmatrixObjectType.tp_init = (initproc)init;
  TmatrixObjectType.tp_new = alloc;

  // finalize creation of the TmatrixObjectType
  PyType_Ready(&TmatrixObjectType);
  // add a TmatrixObjectType to the module
  Py_INCREF(&TmatrixObjectType);
  PyModule_AddObject(module, "TMatrix",
                     reinterpret_cast<PyObject *>(&TmatrixObjectType));
}

#endif // PYTMATRIX_HPP
