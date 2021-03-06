/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TMatrixResource.hpp
 *
 * @brief TMatrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TMATRIXRESOURCE_HPP
#define TMATRIXRESOURCE_HPP

#include "Configuration.hpp"
#include "ConvergedSizeResources.hpp"
#include "GaussBasedResources.hpp"
#include "InteractionResource.hpp"
#include "Matrix.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "WignerDResources.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief TMatrix.
 *
 * Note that, despite the name, the T-matrix does not actually act as a single
 * QuickSched Resource. Instead, it uses its own iternal array of resources
 * to restrict access to the elements corresponding to individual @f$m@f$
 * values. Tasks accessing these elements should call the get_m_resource()
 * member function to obtain the corresponding resource. The @f$m=0@f$ resource
 * additionally also controls access to the internal @f$n_{max}@f$ and
 * @f$n_{GL}@f$ counters and the variables used to check for convergence of
 * the T-matrix calculation procedure.
 */
class TMatrixResource : public Resource {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

  /*! @brief Grant access to @f$m>0@f$ computation task. */
  friend class TMatrixMAllTask;

  /*! @brief Grant access to the scattering and extinction computation task. */
  friend class TMatrixQTask;

  /*! @brief Grant access to the alignment averaging task. */
  friend class AlignmentAverageTask;

  /*! @brief Grant access to the reset task. */
  friend class ResetTMatrixResourceTask;

private:
  /*! @brief Actual size of the matrix stored in this resource. */
  uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points used to compute this
   *  T-matrix. */
  uint_fast32_t _ngauss;

  /*! @brief T-matrix itself. Is in fact a @f$n_{max}+1@f$ element vector for
   *  which every element is a @f$2n_{max}\times{}2n_{max}@f$ matrix. */
  std::vector<Matrix<std::complex<float_type>>> _T;

  /*! @brief Resources for the various @f$m@f$ values stored in the T-matrix. */
  std::vector<Resource> _m_resources;

  /*! @brief Scattering coefficient for the current T-matrix. */
  float_type _Qscattering;

  /*! @brief Exctinction coefficient for the current T-matrix. */
  float_type _Qextinction;

  /*! @brief Relative difference between the scattering coefficient for the
   *  current matrix and that one for the previous T-matrix held by this
   *  resource (if any). */
  float_type _dQscattering;

  /*! @brief Relative difference between the exctinction coefficient for the
   *  current matrix and that one for the previous T-matrix held by this
   *  resource (if any). */
  float_type _dQextinction;

public:
  /**
   * @brief Constructor.
   *
   * Note that we use the fact that the relative differences for the scattering
   * and extinction coefficients are positive to use them as flags: -2 means
   * the TMatrixResource does not contain a valid T-matrix and Q coefficients
   * yet, -1 means this is the first T-matrix held by the resource, so there
   * is no Q coefficient history yet, while positive values correspond to normal
   * operation.
   *
   * @param maximum_order Maximum order of T matrix that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline TMatrixResource(const uint_fast32_t maximum_order)
      : _nmax(0), _ngauss(0), _m_resources(maximum_order + 1), _Qscattering(0.),
        _Qextinction(0.), _dQscattering(-2.), _dQextinction(-2.) {

    _T.reserve(maximum_order + 1);
    for (uint_fast32_t m = 0; m < maximum_order + 1; ++m) {
      const uint_fast32_t nm = maximum_order + 1 - m;
      _T.push_back(Matrix<std::complex<float_type>>(2 * nm, 2 * nm));
    }
  }

  /**
   * @brief Get the size in memory of a hypothetical TMatrixResource object with
   * the given parameters.
   *
   * @param maximum_order Maximum order of the hypothetical object.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t maximum_order) {
    // storage space for class variables
    size_t size = sizeof(TMatrixResource);
    // storage space occupied by T matrix elements
    for (uint_fast32_t m = 0; m < maximum_order + 1; ++m) {
      const uint_fast32_t nm = maximum_order + 1 - m;
      size += 8 * nm * nm * sizeof(float_type);
    }
    return size;
  }

  /**
   * @brief Get the resource corresponding to the given @f$m@f$ value.
   *
   * @param m @f$m@f$ value.
   * @return Corresponding Resource.
   */
  inline Resource &get_m_resource(const uint_fast32_t m) {
    return _m_resources[m];
  }

  /**
   * @brief Get a specific element of the T-matrix.
   *
   * The T-matrix consists of four blocks:
   * @f[
   *    T = \begin{pmatrix}
   *      T^{(11)} & T^{(12)} \\
   *      T^{(21)} & T^{(22)}
   *    \end{pmatrix}.
   * @f]
   * In principle, each of these blocks is an @f$L_{max}\times{}L_{max}@f$
   * matrix indexed using the combined index
   * @f[
   *    l = n (n+1) + m,
   * @f]
   * where @f$n = 1...n_{max}@f$ and @f$m=-n...n@f$. This notation is based on
   * Tsang, Kong & Shin, 1984, Radio Science, 19, 629
   * (https://doi.org/10.1029/RS019i002p00629).
   *
   * This function returns the element @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   *
   * However, since for spheroidal particles @f$T_{nmn'm'} \sim{}
   * \delta{}_{mm'}@f$, we can greatly reduce the storage space of the T matrix
   * by storing it as @f$n_{max}+1@f$ individual matrices (one for each value of
   * @f$m@f$). This is in fact what we do.
   *
   * @param i1 Row index of the desired T-matrix quarter, @f$i_1@f$.
   * @param n1 Order of the row index of the element, @f$n_1@f$.
   * @param m1 Degree of the row index of the element, @f$m_1@f$.
   * @param i2 Column index of the desired T-matrix quarter, @f$i_2@f$.
   * @param n2 Order of the column index of the element, @f$n_2@f$.
   * @param m2 Degree of the column index of the element, @f$m_2@f$.
   * @return Corresponding element, @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   */
  inline const std::complex<float_type> &
  operator()(const uint_fast8_t i1, const uint_fast32_t n1,
             const uint_fast32_t m1, const uint_fast8_t i2,
             const uint_fast32_t n2, const uint_fast32_t m2) const {

    ctm_assert(m1 == m2);
    ctm_assert(m1 <= _nmax);
    ctm_assert(m2 <= _nmax);
    ctm_assert(n1 > 0);
    ctm_assert(n1 <= _nmax);
    ctm_assert(n2 > 0);
    ctm_assert(n2 <= _nmax);
    if (m1 > 0) {
      const uint_fast32_t nm = _nmax + 1 - m1;
      return _T[m1](i1 * nm + n1 - m1, i2 * nm + n2 - m2);
    } else {
      return _T[m1](i1 * _nmax + n1 - m1 - 1, i2 * _nmax + n2 - m2 - 1);
    }
  }

  /**
   * @brief Set the T-matrix element with the given indices to the given value.
   *
   * See operator()() for detailed documentation about how elements are stored.
   *
   * @param i1 Row index of the desired T-matrix quarter, @f$i_1@f$.
   * @param n1 Order of the row index of the element, @f$n_1@f$.
   * @param m1 Degree of the row index of the element, @f$m_1@f$.
   * @param i2 Column index of the desired T-matrix quarter, @f$i_2@f$.
   * @param n2 Order of the column index of the element, @f$n_2@f$.
   * @param m2 Degree of the column index of the element, @f$m_2@f$.
   * @param value New value for the corresponding element,
   * @f$T^{(i_1i_2)}_{n_1n_2m_1m_2}@f$.
   */
  inline void set_element(const uint_fast8_t i1, const uint_fast32_t n1,
                          const uint_fast32_t m1, const uint_fast8_t i2,
                          const uint_fast32_t n2, const uint_fast32_t m2,
                          const std::complex<float_type> &value) {

    ctm_assert(m1 == m2);
    ctm_assert(m1 <= _nmax);
    ctm_assert(m2 <= _nmax);
    ctm_assert(n1 > 0);
    ctm_assert(n1 <= _nmax);
    ctm_assert(n2 > 0);
    ctm_assert(n2 <= _nmax);
    if (m1 > 0) {
      const uint_fast32_t nm = _nmax + 1 - m1;
      _T[m1](i1 * nm + n1 - m1, i2 * nm + n2 - m2) = value;
    } else {
      _T[m1](i1 * _nmax + n1 - m1 - 1, i2 * _nmax + n2 - m2 - 1) = value;
    }
  }

  /**
   * @brief Get the scattering coefficient for the T-matrix.
   *
   * @return Scattering coefficient.
   */
  inline float_type get_scattering_coefficient() const { return _Qscattering; }

  /**
   * @brief Get the extinction coefficient for the T-matrix.
   *
   * @return Extinction coefficient.
   */
  inline float_type get_extinction_coefficient() const { return _Qextinction; }

  /**
   * @brief Get the order of the T-matrix.
   *
   * @return Maximum order, @f$n_{max}@f$.
   */
  inline uint_fast32_t get_nmax() const { return _nmax; }

  /**
   * @brief Get the number of Gauss-Legendre quadrature points for the T-matrix.
   *
   * @return Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline uint_fast32_t get_ngauss() const { return _ngauss; }
};

/**
 * @brief Auxiliary data space for the T matrix calculation.
 */
class TMatrixAuxiliarySpace {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

  /*! @brief Grant access to @f$m>0@f$ computation task. */
  friend class TMatrixMAllTask;

  /*! @brief Grant access to the alignment averaging task. */
  friend class AlignmentAverageTask;

private:
  /*! @brief Q matrix. */
  Matrix<std::complex<float_type>> _Q;

  /*! @brief Regular Q matrix. */
  Matrix<std::complex<float_type>> _RgQ;

  /*! @brief Pivot array for Q matrix inversion. */
  std::vector<uint_fast32_t> _pivot_array;

  /*! @brief Work array for Q matrix inversion. */
  std::vector<std::complex<float_type>> _work;

  /*! @brief Temporary Wigner D function array. */
  std::vector<float_type> _wigner_d;

  /*! @brief Temporary Wigner D function derivative array. */
  std::vector<float_type> _dwigner_d;

  /*! @brief Temporary Clebsch-Gordan coefficient array. */
  std::vector<float_type> _clebsch_gordan;

public:
  /**
   * @brief Constructor.
   *
   * @param maximum_order Maximum order of T matrix that will ever be requested
   * from this resource, used to initialise the internal storage space.
   */
  inline TMatrixAuxiliarySpace(const uint_fast32_t maximum_order)
      : _Q(2 * maximum_order, 2 * maximum_order),
        _RgQ(2 * maximum_order, 2 * maximum_order),
        _pivot_array(2 * maximum_order), _work(2 * maximum_order),
        _wigner_d(maximum_order), _dwigner_d(maximum_order),
        _clebsch_gordan((maximum_order + 1) * (2 * maximum_order + 1)) {}

  /**
   * @brief Get the size in memory of a hypothetical TMatrixAuxiliarySpace
   * object with the given parameters.
   *
   * @param maximum_order Maximum order of the hypothetical object.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t maximum_order) {
    // storage space for class variables
    size_t size = sizeof(TMatrixAuxiliarySpace);
    // Q
    size += 8 * maximum_order * maximum_order * sizeof(float_type);
    // RgQ
    size += 8 * maximum_order * maximum_order * sizeof(float_type);
    // pivot
    size += 2 * maximum_order * sizeof(uint_fast32_t);
    // work
    size += 4 * maximum_order * sizeof(float_type);
    // wigner_d
    size += 2 * maximum_order * sizeof(float_type);
    // clesch_gordan
    size += (maximum_order + 1) * (2 * maximum_order + 1) * sizeof(float_type);
    return size;
  }

  /**
   * @brief Clear the entire contents of the matrices.
   */
  inline void reset() {
    _Q.reset();
    _RgQ.reset();
    for (uint_fast32_t i = 0; i < _pivot_array.size(); ++i) {
      _pivot_array[i] = 0;
      _work[i] = std::complex<float_type>(0.);
    }
    for (uint_fast32_t i = 0; i < _wigner_d.size(); ++i) {
      _wigner_d[i] = 0.;
      _dwigner_d[i] = 0.;
    }
    for (uint_fast32_t i = 0; i < _clebsch_gordan.size(); ++i) {
      _clebsch_gordan[i] = 0.;
    }
  }
};

/**
 * @brief Manager for TMatrixAuxiliary space classes.
 */
class TMatrixAuxiliarySpaceManager {
private:
  /*! @brief Vector of preallocated auxiliary spaces. */
  std::vector<TMatrixAuxiliarySpace *> _spaces;

public:
  /**
   * @brief Constructor.
   *
   * @param num_threads Number of threads that is used (and number of
   * auxiliary spaces to allocate).
   * @param maximum_order Maximum order of any T-matrix that will be computed,
   * used to set the size of the auxiliary space.
   */
  inline TMatrixAuxiliarySpaceManager(const int_fast32_t num_threads,
                                      const uint_fast32_t maximum_order)
      : _spaces(num_threads, nullptr) {
    for (int_fast32_t i = 0; i < num_threads; ++i) {
      _spaces[i] = new TMatrixAuxiliarySpace(maximum_order);
    }
  }

  /**
   * @brief Destructor.
   *
   * Deallocate the auxiliary spaces.
   */
  inline ~TMatrixAuxiliarySpaceManager() {
    for (uint_fast32_t i = 0; i < _spaces.size(); ++i) {
      delete _spaces[i];
    }
  }

  /**
   * @brief Get the size in memory of a hypothetical
   * TMatrixAuxiliarySpaceManager object with the given parameters.
   *
   * @param num_threads Number of threads or number of auxiliary spaces.
   * @param maximum_order Maximum order of the hypothetical spaces.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const int_fast32_t num_threads,
                                       const uint_fast32_t maximum_order) {
    // storage space for class variables
    size_t size = sizeof(TMatrixAuxiliarySpaceManager);
    size += num_threads * TMatrixAuxiliarySpace::get_memory_size(maximum_order);
    return size;
  }

  /**
   * @brief Get a reference to the auxiliary space for the thread with the
   * given ID.
   *
   * @param thread_id Thread ID.
   * @return Corresponding auxiliary space.
   */
  inline TMatrixAuxiliarySpace &get_space(int_fast32_t thread_id) {
    return *_spaces[thread_id];
  }
};

/**
 * @brief Compute the @f$m=0@f$ part of the T-matrix for a given @f$n_{max}@f$.
 */
class TMatrixM0Task : public Task {
private:
  /*! @brief Desired accuracy for the T-matrix calculation. The task will only
   *  recompute the T-matrix elements if the Q coefficients in the
   *  TMatrixResource have not reached this accuracy yet. */
  const float_type _tolerance;

  /*! @brief Maximum order of the matrix, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  const uint_fast32_t _ngauss;

  /*! @brief Precomputed @f$n@f$ factors (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Precomputed Gauss-Legendre quadrature points (read only). */
  const GaussBasedResources &_quadrature_points;

  /*! @brief Precomputed geometry specific quadrature points (read only). */
  const ParticleGeometryResource &_geometry;

  /*! @brief Interaction variables (read only). */
  const InteractionVariables _interaction_variables;

  /*! @brief Precomputed interaction specific quadrature points (read only). */
  const InteractionResource &_interaction;

  /*! @brief Precomputed @f$m=0@f$ Wigner D functions (read only). */
  const WignerDResources &_wigner;

  /*! @brief Auxiliary space manager used to obtain space to store intermediate
   *  calculations. */
  TMatrixAuxiliarySpaceManager &_aux_manager;

  /*! @brief TMatrix space containing the result. */
  TMatrixResource &_Tmatrix;

  /*! @brief Resource containing the converged order and number of quadrature
   *  points. */
  ConvergedSizeResources &_converged_size;

  /*! @brief Resource that guarantees unique access to the @f$m=0@f$ T-matrix
   *  elements. */
  const Resource &_m_resource;

public:
  /**
   * @brief Constructor.
   *
   * @param tolerance Desired accuracy for the T-matrix calculation. The task
   * will only recompute the T-matrix elements if the Q coefficients in the
   * TMatrixResource have not reached this accuracy yet.
   * @param nmax Maximum order of the matrix, @f$n_{max}@f$.
   * @param ngauss Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   * @param nfactors Precomputed @f$n@f$ factors (read only).
   * @param quadrature_points Precomputed Gauss-Legendre quadrature points (read
   * only).
   * @param geometry Precomputed geometry specific quadrature points (read
   * only).
   * @param interaction_variables Interaction variables (read only).
   * @param interaction Precomputed interaction specific quadrature points (read
   * only).
   * @param wigner Precomputed @f$m=0@f$ Wigner D functions (read only).
   * @param aux_manager Auxiliary space manager used to obtain space to store
   * intermediate calculations.
   * @param Tmatrix TMatrix space containing the result.
   * @param converged_size Resource containing the converged order and number of
   * quadrature points.
   * @param m_resource Resource that guarantees unique access to the @f$m=0@f$
   * values of the T-matrix.
   */
  inline TMatrixM0Task(
      const float_type tolerance, const uint_fast32_t nmax,
      const uint_fast32_t ngauss, const NBasedResources &nfactors,
      const GaussBasedResources &quadrature_points,
      const ParticleGeometryResource &geometry,
      const InteractionVariables &interaction_variables,
      const InteractionResource &interaction, const WignerDResources &wigner,
      TMatrixAuxiliarySpaceManager &aux_manager, TMatrixResource &Tmatrix,
      ConvergedSizeResources &converged_size, const Resource &m_resource)
      : _tolerance(tolerance), _nmax(nmax), _ngauss(ngauss),
        _nfactors(nfactors), _quadrature_points(quadrature_points),
        _geometry(geometry), _interaction_variables(interaction_variables),
        _interaction(interaction), _wigner(wigner), _aux_manager(aux_manager),
        _Tmatrix(Tmatrix), _converged_size(converged_size),
        _m_resource(m_resource) {}

  virtual ~TMatrixM0Task() {}

  /**
   * @brief Get the size in memory of a hypothetical TMatrixM0Task object with
   * the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(TMatrixM0Task);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _m_resource, true);
    quicksched.link_task_and_resource(*this, _converged_size, true);

    // read access
    quicksched.link_task_and_resource(*this, _nfactors, false);
    quicksched.link_task_and_resource(*this, _quadrature_points, false);
    quicksched.link_task_and_resource(*this, _geometry, false);
    quicksched.link_task_and_resource(*this, _interaction_variables, false);
    quicksched.link_task_and_resource(*this, _interaction, false);
    quicksched.link_task_and_resource(*this, _wigner, false);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 2.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 2; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 6.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 6; }

  /**
   * @brief Compute the @f$m=0@f$ elements of the T-matrix.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    // check that the T-matrix resource was actually reset properly
    ctm_assert(_converged_size._iteration > 0 || _Tmatrix._nmax == 0);

    // check if we need to do something
    if (_converged_size.is_converged()) {
      // tolerance was already reached; abort task
      return;
    }

    TMatrixAuxiliarySpace &aux = _aux_manager.get_space(thread_id);

    aux.reset();

    std::vector<float_type> &wigner_d = aux._wigner_d;
    std::vector<float_type> &dwigner_d = aux._dwigner_d;

    const std::complex<float_type> k2mr =
        _interaction_variables.get_material_wavenumber_times_wavenumber();
    const float_type k2 = _interaction_variables.get_wavenumber_squared();

    for (uint_fast32_t ig = 1; ig < _ngauss + 1; ++ig) {

      _wigner.get_wigner_functions(0, ig - 1, wigner_d, dwigner_d);

      for (uint_fast32_t n1 = 1; n1 < _nmax + 1; ++n1) {

        // n1 * (n1 + 1)
        const float_type n1n1p1 = n1 * (n1 + 1.);

        const uint_fast32_t i1 = n1;
        const uint_fast32_t ii1 = i1 + _nmax;

        for (uint_fast32_t n2 = 1; n2 < _nmax + 1; ++n2) {

          // n2 * (n2 + 1)
          const float_type n2n2p1 = n2 * (n2 + 1.);

          const uint_fast32_t i2 = n2;
          const uint_fast32_t ii2 = i2 + _nmax;

          // filter out half the components because of symmetry
          if ((n1 + n2) % 2 == 0) {
            const float_type wigner_n1 = wigner_d[n1 - 1];
            const float_type dwigner_n1 = dwigner_d[n1 - 1];
            const float_type wigner_n2 = wigner_d[n2 - 1];
            const float_type dwigner_n2 = dwigner_d[n2 - 1];

            const float_type wn1dwn2 = wigner_n1 * dwigner_n2;
            const float_type dwn1wn2 = dwigner_n1 * wigner_n2;
            const float_type dwn1dwn2 = dwigner_n1 * dwigner_n2;

            const float_type jkrn1 = _interaction.get_jkr(ig - 1, n1);
            const float_type ykrn1 = _interaction.get_ykr(ig - 1, n1);
            // spherical Hankel function of the first kind
            const std::complex<float_type> hkrn1(jkrn1, ykrn1);
            const float_type djkrn1 = _interaction.get_djkr(ig - 1, n1);
            const float_type dykrn1 = _interaction.get_dykr(ig - 1, n1);
            const std::complex<float_type> dhkrn1(djkrn1, dykrn1);
            const std::complex<float_type> jkrmrn2 =
                _interaction.get_jkrmr(ig - 1, n2);
            const std::complex<float_type> djkrmrn2 =
                _interaction.get_djkrmr(ig - 1, n2);

            const std::complex<float_type> c1 = jkrmrn2 * jkrn1;
            const std::complex<float_type> b1 = jkrmrn2 * hkrn1;

            const std::complex<float_type> c2 = jkrmrn2 * djkrn1;
            const std::complex<float_type> b2 = jkrmrn2 * dhkrn1;

            const float_type krinvi = _interaction.get_krinv(ig - 1);
            const std::complex<float_type> c3 = krinvi * c1;
            const std::complex<float_type> b3 = krinvi * b1;

            const std::complex<float_type> c4 = jkrn1 * djkrmrn2;
            const std::complex<float_type> b4 = hkrn1 * djkrmrn2;

            const std::complex<float_type> krmrinvi =
                _interaction.get_krmrinv(ig - 1);
            const std::complex<float_type> c5 = c1 * krmrinvi;
            const std::complex<float_type> b5 = b1 * krmrinvi;

            // note that r2 here actually is (r2/R_V^2)
            // this is okay, since this factor appears in all elements of the
            // Q matrix, so it will divide out when multiplying with the inverse
            const float_type wr2i = _quadrature_points.get_weight(ig - 1) *
                                    _geometry.get_r2(ig - 1);
            const float_type dr_over_ri = _geometry.get_dr_over_r(ig - 1);

            // prefactor sqrt{(2n1+1)*(2n2+1)/[n1*(n1+1)*n2*(n2+1)]}
            const float_type an12 = _nfactors.get_ann(n1, n2);

            const float_type f1 = wr2i * dwn1dwn2;
            const float_type f2 = wr2i * dr_over_ri * n1n1p1 * wn1dwn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*hn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * hn1(kr) / kr
            const std::complex<float_type> this_J12 =
                -an12 * (f1 * b2 + f2 * b3);
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn2(krmr) *
            //    ([kr*jn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * jn1(kr) / kr
            const std::complex<float_type> this_RgJ12 =
                -an12 * (f1 * c2 + f2 * c3);

            const float_type f3 = wr2i * dr_over_ri * n2n2p1 * dwn1wn2;
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * hn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * hn1(kr) / krmr
            const std::complex<float_type> this_J21 =
                an12 * (f1 * b4 + f3 * b5);
            // r^2 * ddn1_0m/dtheta * ddn2_0m/dtheta * jn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * jn1(kr) / krmr
            const std::complex<float_type> this_RgJ21 =
                an12 * (f1 * c4 + f3 * c5);

            aux._Q(i1 - 1, i2 - 1) += k2mr * this_J21 + k2 * this_J12;
            aux._RgQ(i1 - 1, i2 - 1) += k2mr * this_RgJ21 + k2 * this_RgJ12;

            aux._Q(ii1 - 1, ii2 - 1) += k2mr * this_J12 + k2 * this_J21;
            aux._RgQ(ii1 - 1, ii2 - 1) += k2mr * this_RgJ12 + k2 * this_RgJ21;
          }
        }
      }
    }

    // func_TT
    const uint_fast32_t nmax2 = 2 * _nmax;
    aux._Q.plu_inverse(nmax2, &aux._pivot_array[0], aux._pivot_array.size(),
                       &aux._work[0], aux._work.size());

    for (uint_fast32_t i = 0; i < _nmax; ++i) {
      for (uint_fast32_t j = 0; j < _nmax; ++j) {
        _Tmatrix._T[0](i, j) = 0.;
        _Tmatrix._T[0](_nmax + i, j) = 0.;
        _Tmatrix._T[0](i, _nmax + j) = 0.;
        _Tmatrix._T[0](_nmax + i, _nmax + j) = 0.;
        for (uint_fast32_t k = 0; k < nmax2; ++k) {
          _Tmatrix._T[0](i, j) -= aux._RgQ(i, k) * aux._Q(k, j);
          _Tmatrix._T[0](_nmax + i, j) -= aux._RgQ(_nmax + i, k) * aux._Q(k, j);
          _Tmatrix._T[0](i, _nmax + j) -= aux._RgQ(i, k) * aux._Q(k, _nmax + j);
          _Tmatrix._T[0](_nmax + i, _nmax + j) -=
              aux._RgQ(_nmax + i, k) * aux._Q(k, _nmax + j);
        }
      }
    }

    // set the order and number of quadrature points of the T-matrix
    _Tmatrix._nmax = _nmax;
    _Tmatrix._ngauss = _ngauss;

    _converged_size._nmax = _nmax;
    _converged_size._ngauss = _ngauss;
    _converged_size._quadrature_points = &_quadrature_points;
    _converged_size._geometry = &_geometry;
    _converged_size._wigner_d = &_wigner;
    ++_converged_size._iteration;

    // update Q coefficients
    const float_type old_Qscattering = _Tmatrix._Qscattering;
    const float_type old_Qextinction = _Tmatrix._Qextinction;
    _Tmatrix._Qscattering = 0.;
    _Tmatrix._Qextinction = 0.;
    for (uint_fast32_t n = 1; n < _nmax + 1; ++n) {
      const float_type dn1 = 2. * n + 1.;
      _Tmatrix._Qscattering +=
          dn1 * (norm(_Tmatrix._T[0](n - 1, n - 1)) +
                 norm(_Tmatrix._T[0](_nmax + n - 1, _nmax + n - 1)));
      _Tmatrix._Qextinction +=
          dn1 * (_Tmatrix._T[0](n - 1, n - 1).real() +
                 _Tmatrix._T[0](_nmax + n - 1, _nmax + n - 1).real());
    }

    if (_Tmatrix._dQscattering == -2.) {
      _Tmatrix._dQscattering = -1.;
      _Tmatrix._dQextinction = -1.;
    } else {
      _Tmatrix._dQscattering = fabs((old_Qscattering - _Tmatrix._Qscattering) /
                                    _Tmatrix._Qscattering);
      _Tmatrix._dQextinction = fabs((old_Qextinction - _Tmatrix._Qextinction) /
                                    _Tmatrix._Qextinction);
    }

    if (_Tmatrix._dQscattering > 0.) {
      _converged_size._dQs = _Tmatrix._dQscattering;
      _converged_size._dQe = _Tmatrix._dQextinction;
      if (_Tmatrix._dQscattering <= _tolerance &&
          _Tmatrix._dQextinction <= _tolerance) {
        _converged_size._is_converged = true;
      }
    }
  }

  /**
   * @brief Get the computational cost of this task.
   *
   * @return Computational cost.
   */
  virtual int_fast32_t get_cost() const {
    return 93468 * _ngauss * _ngauss - 4358242 * _ngauss + 66752438;
  }
};

/**
 * @brief Compute the @f$m>0@f$ part of the T-matrix for a given @f$n_{max}@f$.
 */
class TMatrixMAllTask : public Task {
private:
  /*! @brief @f$m@f$ value. */
  const uint_fast32_t _m;

  /*! @brief Precomputed @f$n@f$ factors (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Converged T-matrix variables. */
  const ConvergedSizeResources &_converged_size;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief Interaction resources. */
  const InteractionResource &_interaction;

  /*! @brief Auxiliary space manager used to obtain space to store intermediate
   *  calculations. */
  TMatrixAuxiliarySpaceManager &_aux_manager;

  /*! @brief TMatrix space containing the result. */
  TMatrixResource &_Tmatrix;

  /*! @brief Resource that guarantees unique access to the @f$m>0@f$ T-matrix
   *  elements. */
  const Resource &_m_resource;

public:
  /**
   * @brief Constructor.
   *
   * @param m @f$m@f$ value.
   * @param nfactors Precomputed @f$n@f$ factors (read only).
   * @param converged_size Converged T-matrix variables (read only).
   * @param interaction_variables Interaction variables (read only).
   * @param interaction Interaction resources (read only).
   * @param aux_manager Auxiliary space manager used to obtain space to store
   * intermediate calculations.
   * @param Tmatrix TMatrix space containing the result.
   * @param m_resource Resource that guarantees unique access to the @f$m=0@f$
   * values of the T-matrix.
   */
  inline TMatrixMAllTask(const uint_fast32_t m, const NBasedResources &nfactors,
                         const ConvergedSizeResources &converged_size,
                         const InteractionVariables &interaction_variables,
                         const InteractionResource &interaction,
                         TMatrixAuxiliarySpaceManager &aux_manager,
                         TMatrixResource &Tmatrix, const Resource &m_resource)
      : _m(m), _nfactors(nfactors), _converged_size(converged_size),
        _interaction_variables(interaction_variables),
        _interaction(interaction), _aux_manager(aux_manager), _Tmatrix(Tmatrix),
        _m_resource(m_resource) {}

  virtual ~TMatrixMAllTask() {}

  /**
   * @brief Get the size in memory of a hypothetical TMatrixMAllTask object with
   * the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(TMatrixMAllTask);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _m_resource, true);

    // read access
    quicksched.link_task_and_resource(*this, _nfactors, false);
    quicksched.link_task_and_resource(*this, _converged_size, false);
    quicksched.link_task_and_resource(*this, _interaction_variables, false);
    quicksched.link_task_and_resource(*this, _interaction, false);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 4.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 4; }

  /**
   * @brief Compute the @f$m>0@f$ elements of the T-matrix.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    // make sure the T-matrix was actually converged
    ctm_assert_message(
        _converged_size.is_converged(),
        "size: %g, wavelength: %g, axis_ratio: %g, dQs: %g, dQe: %g",
        double(_interaction_variables.get_equal_volume_radius()),
        double(2. * M_PI / _interaction_variables.get_wavenumber()),
        double(_converged_size.get_geometry()->get_axis_ratio()),
        double(_converged_size.get_dQs()), double(_converged_size.get_dQe()));

    // since we don't know how many elements the T matrix will have before
    // we create the tasks, we might have tasks for elements of the T matrix
    // that we don't need. We make sure they don't waste computing time.
    const uint_fast32_t nmax = _converged_size.get_nmax();
    if (_m > nmax) {
      return;
    }

    TMatrixAuxiliarySpace &aux = _aux_manager.get_space(thread_id);

    aux.reset();

    const uint_fast32_t ngauss = _converged_size.get_ngauss();
    const GaussBasedResources &quadrature_points =
        *_converged_size.get_quadrature_points();
    const ParticleGeometryResource &geometry = *_converged_size.get_geometry();
    const WignerDResources &wigner_d = *_converged_size.get_wigner_d();

    const float_type m2 = _m * _m;
    const uint_fast32_t nm = nmax + 1 - _m;
    const uint_fast32_t nm2 = 2 * nm;

    std::vector<float_type> &wigner_ds = aux._wigner_d;
    std::vector<float_type> &dwigner_ds = aux._dwigner_d;

    const std::complex<float_type> k2mr =
        _interaction_variables.get_material_wavenumber_times_wavenumber();
    const float_type k2 = _interaction_variables.get_wavenumber_squared();
    const std::complex<float_type> icompl(0., 1.);

    for (uint_fast32_t ig = 0; ig < ngauss; ++ig) {
      wigner_d.get_wigner_functions(_m, ig, wigner_ds, dwigner_ds);
      for (uint_fast32_t n1 = _m; n1 < nmax + 1; ++n1) {

        // n1 * (n1 + 1)
        const float_type n1n1p1 = n1 * (n1 + 1.);

        const uint_fast32_t i1 = n1 + 1 - _m;
        const uint_fast32_t ii1 = i1 + nm;

        for (uint_fast32_t n2 = _m; n2 < nmax + 1; ++n2) {

          // n2 * (n2 + 1)
          const float_type n2n2p1 = n2 * (n2 + 1.);

          const uint_fast32_t i2 = n2 + 1 - _m;
          const uint_fast32_t ii2 = i2 + nm;

          const int_fast8_t sign = ((n1 + n2) % 2 == 0) ? 1 : -1;
          const float_type wigner_n1 = wigner_ds[n1 - 1];
          const float_type dwigner_n1 = dwigner_ds[n1 - 1];
          const float_type wigner_n2 = wigner_ds[n2 - 1];
          const float_type dwigner_n2 = dwigner_ds[n2 - 1];

          const float_type wn1wn2 = wigner_n1 * wigner_n2;
          const float_type wn1dwn2 = wigner_n1 * dwigner_n2;
          const float_type dwn1wn2 = dwigner_n1 * wigner_n2;

          const float_type jkrn1 = _interaction.get_jkr(ig, n1);
          const float_type ykrn1 = _interaction.get_ykr(ig, n1);
          // spherical Hankel function of the first kind
          const std::complex<float_type> hkrn1(jkrn1, ykrn1);
          const float_type djkrn1 = _interaction.get_djkr(ig, n1);
          const float_type dykrn1 = _interaction.get_dykr(ig, n1);
          const std::complex<float_type> dhkrn1(djkrn1, dykrn1);
          const std::complex<float_type> jkrmrn2 =
              _interaction.get_jkrmr(ig, n2);
          const std::complex<float_type> djkrmrn2 =
              _interaction.get_djkrmr(ig, n2);

          const std::complex<float_type> c1 = jkrmrn2 * jkrn1;
          const std::complex<float_type> b1 = jkrmrn2 * hkrn1;

          const std::complex<float_type> c2 = jkrmrn2 * djkrn1;
          const std::complex<float_type> b2 = jkrmrn2 * dhkrn1;

          const float_type krinvi = _interaction.get_krinv(ig);

          const std::complex<float_type> c4 = djkrmrn2 * jkrn1;
          const std::complex<float_type> b4 = djkrmrn2 * hkrn1;

          const std::complex<float_type> krmrinvi =
              _interaction.get_krmrinv(ig);

          const float_type dr_over_ri = geometry.get_dr_over_r(ig);

          // prefactor sqrt{(2n1+1)*(2n2+1)/[n1*(n1+1)*n2*(n2+1)]}
          const float_type an12 = _nfactors.get_ann(n1, n2);

          if (sign < 0) {
            // note that r2 here actually is (r2/R_V^2)
            // this is okay, since this factor appears in all elements of the
            // Q matrix, so it will divide out when multiplying with the inverse
            const float_type dsi = quadrature_points.get_sinthetainv(ig) * _m *
                                   quadrature_points.get_weight(ig) *
                                   geometry.get_r2(ig);

            const std::complex<float_type> c6 = djkrmrn2 * djkrn1;
            const std::complex<float_type> b6 = djkrmrn2 * dhkrn1;

            const std::complex<float_type> c7 = c4 * krinvi;
            const std::complex<float_type> b7 = b4 * krinvi;

            const std::complex<float_type> c8 = c2 * krmrinvi;
            const std::complex<float_type> b8 = b2 * krmrinvi;

            const float_type e1 = dsi * (wn1dwn2 + dwn1wn2);
            // (m / sintheta) * jn2(krmr) * hn1(kr) *
            // (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
            const std::complex<float_type> this_J11 = icompl * an12 * (e1 * b1);
            // (m / sintheta) * jn2(krmr) * jn1(kr) *
            // (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
            const std::complex<float_type> this_RgJ11 =
                icompl * an12 * (e1 * c1);

            const float_type factor = dsi * dr_over_ri * wn1wn2;
            const float_type e2 = factor * n1n1p1;
            const float_type e3 = factor * n2n2p1;
            // (m / sintheta) * ([krmr*kn2(krmr)]'/krmr) * ([kr*hn1(kr)]'/kr)
            //    * (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
            // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n1 * (n1 + 1)
            //    * ([krmr*jn2(krmr)]'/krmr) * hn1(kr) / kr
            // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n2 * (n2 + 1)
            //    * jn2(krmr) * ([kr*hn1(kr]'/kr) / krmr
            const std::complex<float_type> this_J22 =
                icompl * an12 * (e1 * b6 + e2 * b7 + e3 * b8);
            // (m / sintheta) * ([krmr*kn2(krmr)]'/krmr) * ([kr*jn1(kr)]'/kr)
            //    * (dn1_0m * ddn2_0m/dtheta + ddn1_0m/dtheta * dn2_0m)
            // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n1 * (n1 + 1)
            //    * ([krmr*jn2(krmr)]'/krmr) * jn1(kr) / kr
            // + (m / sintheta) * dr/rdtheta * dn1_0m * dn2_0m * n2 * (n2 + 1)
            //    * jn2(krmr) * ([kr*jn1(kr]'/kr) / krmr
            const std::complex<float_type> this_RgJ22 =
                icompl * an12 * (e1 * c6 + e2 * c7 + e3 * c8);

            aux._Q(i1 - 1, ii2 - 1) += k2mr * this_J11 + k2 * this_J22;
            aux._RgQ(i1 - 1, ii2 - 1) += k2mr * this_RgJ11 + k2 * this_RgJ22;

            aux._Q(ii1 - 1, i2 - 1) += k2mr * this_J22 + k2 * this_J11;
            aux._RgQ(ii1 - 1, i2 - 1) += k2mr * this_RgJ22 + k2 * this_RgJ11;

          } else {
            // note that r2 here actually is (r2/R_V^2)
            // this is okay, since this factor appears in all elements of the
            // Q matrix, so it will divide out when multiplying with the inverse
            const float_type wr2i =
                quadrature_points.get_weight(ig) * geometry.get_r2(ig);
            const float_type dssi = quadrature_points.get_sintheta2inv(ig) * m2;

            const std::complex<float_type> c3 = krinvi * c1;
            const std::complex<float_type> b3 = krinvi * b1;
            const std::complex<float_type> c5 = c1 * krmrinvi;
            const std::complex<float_type> b5 = b1 * krmrinvi;

            const float_type f1 =
                wr2i * (wn1wn2 * dssi + dwigner_n1 * dwigner_n2);
            const float_type f2 = wr2i * dr_over_ri * n1n1p1 * wn1dwn2;
            // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
            //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn2(krmr) *
            //    ([kr*hn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * hn1(kr) / kr
            const std::complex<float_type> this_J12 =
                -an12 * (f1 * b2 + f2 * b3);
            // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
            //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn2(krmr) *
            //    ([kr*jn1(kr)]'/kr)
            // + r^2 * dr/rdtheta * n1 * (n1 + 1) * dn1_0m * ddn1_0m/dtheta *
            //    jn2(krmr) * jn1(kr) / kr
            const std::complex<float_type> this_RgJ12 =
                -an12 * (f1 * c2 + f2 * c3);

            const float_type f3 = wr2i * dr_over_ri * n2n2p1 * dwn1wn2;
            // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
            //        ddn1_0m/dtheta * ddn2_0m/dtheta) * hn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * hn1(kr) / krmr
            const std::complex<float_type> this_J21 =
                an12 * (f1 * b4 + f3 * b5);
            // r^2 * (m^2 / sintheta * dn1_0m * dn2_0m +
            //        ddn1_0m/dtheta * ddn2_0m/dtheta) * jn1(kr) *
            //    ([krmr*jn2(krmr)]'/krmr)
            // + r^2 * dr/rdtheta * n2 * (n2 + 1) * ddn1_0m/dtheta * dn2_0m *
            //    jn2(krmr) * jn1(kr) / krmr
            const std::complex<float_type> this_RgJ21 =
                an12 * (f1 * c4 + f3 * c5);

            aux._Q(i1 - 1, i2 - 1) += k2mr * this_J21 + k2 * this_J12;
            aux._RgQ(i1 - 1, i2 - 1) += k2mr * this_RgJ21 + k2 * this_RgJ12;

            aux._Q(ii1 - 1, ii2 - 1) += k2mr * this_J12 + k2 * this_J21;
            aux._RgQ(ii1 - 1, ii2 - 1) += k2mr * this_RgJ12 + k2 * this_RgJ21;
          }
        }
      }
    }

    aux._Q.plu_inverse(nm2, &aux._pivot_array[0], aux._pivot_array.size(),
                       &aux._work[0], aux._work.size());

    for (uint_fast32_t i = 0; i < nm; ++i) {
      for (uint_fast32_t j = 0; j < nm; ++j) {
        _Tmatrix._T[_m](i, j) = 0.;
        _Tmatrix._T[_m](nm + i, j) = 0.;
        _Tmatrix._T[_m](i, nm + j) = 0.;
        _Tmatrix._T[_m](nm + i, nm + j) = 0.;
        for (uint_fast32_t k = 0; k < nm2; ++k) {
          _Tmatrix._T[_m](i, j) -= aux._RgQ(i, k) * aux._Q(k, j);
          _Tmatrix._T[_m](nm + i, j) -= aux._RgQ(nm + i, k) * aux._Q(k, j);
          _Tmatrix._T[_m](i, nm + j) -= aux._RgQ(i, k) * aux._Q(k, nm + j);
          _Tmatrix._T[_m](nm + i, nm + j) -=
              aux._RgQ(nm + i, k) * aux._Q(k, nm + j);
        }
      }
    }
  }

  /**
   * @brief Get the computational cost of this task.
   *
   * @return Computational cost.
   */
  virtual int_fast32_t get_cost() const {
    return static_cast<int_fast32_t>(std::round(2.5e8 * std::exp(-0.117 * _m)));
  }
};

/**
 * @brief Task that computes the scattering and extinction coefficients for a
 * complete T-matrix.
 */
class TMatrixQTask : public Task {
private:
  /*! @brief TMatrix space containing the result. */
  TMatrixResource &_Tmatrix;

public:
  /**
   * @brief Constructor.
   *
   * @param Tmatrix TMatrix resource containing the full T-matrix.
   */
  inline TMatrixQTask(TMatrixResource &Tmatrix) : _Tmatrix(Tmatrix) {}

  virtual ~TMatrixQTask() {}

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _Tmatrix, true);
  }

  /**
   * @brief Compute the scattering and extinction coefficients.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    _Tmatrix._Qscattering = 0.;
    for (uint_fast32_t n1 = 1; n1 < _Tmatrix._nmax + 1; ++n1) {
      for (uint_fast32_t n2 = 1; n2 < _Tmatrix._nmax + 1; ++n2) {
        for (uint_fast32_t m1 = 0; m1 < n1 + 1; ++m1) {
          for (uint_fast32_t m2 = 0; m2 < n2 + 1; ++m2) {
            if (m1 == m2) {
              float_type factor;
              if (m1 > 0) {
                factor = 2.;
              } else {
                factor = 1.;
              }
              _Tmatrix._Qscattering +=
                  factor * std::norm(_Tmatrix(0, n1, m1, 0, n2, m2));
              _Tmatrix._Qscattering +=
                  factor * std::norm(_Tmatrix(0, n1, m1, 1, n2, m2));
              _Tmatrix._Qscattering +=
                  factor * std::norm(_Tmatrix(1, n1, m1, 0, n2, m2));
              _Tmatrix._Qscattering +=
                  factor * std::norm(_Tmatrix(1, n1, m1, 1, n2, m2));
            }
          }
        }
      }
    }

    _Tmatrix._Qextinction = 0.;
    for (uint_fast32_t n1 = 1; n1 < _Tmatrix._nmax + 1; ++n1) {
      for (uint_fast32_t m1 = 0; m1 < n1 + 1; ++m1) {
        float_type factor;
        if (m1 > 0) {
          factor = 2.;
        } else {
          factor = 1.;
        }
        _Tmatrix._Qextinction += factor * _Tmatrix(0, n1, m1, 0, n1, m1).real();
        _Tmatrix._Qextinction += factor * _Tmatrix(0, n1, m1, 1, n1, m1).real();
        _Tmatrix._Qextinction += factor * _Tmatrix(1, n1, m1, 0, n1, m1).real();
        _Tmatrix._Qextinction += factor * _Tmatrix(1, n1, m1, 1, n1, m1).real();
      }
    }
  }
};

/**
 * @brief Reset a TMatrixResource so that it is ready to be reused.
 */
class ResetTMatrixResourceTask : public Task {
private:
  /*! @brief TMatrixResource to reset. */
  TMatrixResource &_Tmatrix;

public:
  /**
   * @brief Constructor.
   *
   * @param Tmatrix TMatrixResource to reset.
   */
  inline ResetTMatrixResourceTask(TMatrixResource &Tmatrix)
      : _Tmatrix(Tmatrix) {}

  virtual ~ResetTMatrixResourceTask() {}

  /**
   * @brief Get the size in memory of a hypothetical ResetTMatrixResourceTask
   * object with the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(ResetTMatrixResourceTask);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _Tmatrix, true);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 0.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 0; }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    _Tmatrix._nmax = 0;
    _Tmatrix._ngauss = 0;
    _Tmatrix._Qscattering = 0.;
    _Tmatrix._Qextinction = 0.;
    _Tmatrix._dQscattering = -2.;
    _Tmatrix._dQextinction = -2.;
    for (uint_fast32_t m = 0; m < _Tmatrix._m_resources.size(); ++m) {
      _Tmatrix._T[m].reset();
    }
  }
};

#endif // TMATRIXRESOURCE_HPP
