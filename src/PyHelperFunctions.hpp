/**
 * @file PyHelperFunctions.hpp
 *
 * @brief Header library containing helper functions to deal with the Python C
 * API.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PYHELPERFUNCTIONS_HPP
#define PYHELPERFUNCTIONS_HPP

#include "Configuration.hpp"
#include "Error.hpp"

#include <Python.h>
#include <vector>
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unpack the given 0D or 1D NumPy array into a std::vector with the
 * same number of elements.
 *
 * @param numpy_array Input NumPy array.
 * @return Output std::vector.
 * @tparam DATA_TYPE Data type of the elements in the result vector.
 * @tparam INPUT_TYPE Data type to use when reading values from the NumPy array.
 */
template <typename DATA_TYPE, typename INPUT_TYPE = DATA_TYPE>
inline std::vector<DATA_TYPE> unpack_numpy_array(PyArrayObject *numpy_array) {

  // determine the size of the array
  npy_intp array_size;
  const npy_intp array_ndim = PyArray_NDIM(numpy_array);
  if (array_ndim > 1) {
    ctm_warning("Wrong shape for input array!");
    return std::vector<DATA_TYPE>();
  }
  if (array_ndim > 0) {
    const npy_intp *array_dims = PyArray_DIMS(numpy_array);
    array_size = array_dims[0];
  } else {
    array_size = 1;
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    numpy_array = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(numpy_array, &newdimsobj, NPY_ANYORDER));
  }

  std::vector<DATA_TYPE> array(array_size);
  for (npy_intp i = 0; i < array_size; ++i) {
    array[i] =
        *(reinterpret_cast<INPUT_TYPE *>(PyArray_GETPTR1(numpy_array, i)));
  }

  return array;
}

#endif // PYHELPERFUNCTIONS_HPP
