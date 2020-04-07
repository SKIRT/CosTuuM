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
 * @file ConvergedSizeResources.hpp
 *
 * @brief Maximum order and number of Gauss-Legendre quadrature points for a
 * converged T-matrix.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CONVERGEDSIZERESOURCES_HPP
#define CONVERGEDSIZERESOURCES_HPP

#include "Configuration.hpp"
#include "Matrix.hpp"
#include "QuickSchedWrapper.hpp"

#include <cmath>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

class GaussBasedResources;
class InteractionResource;
class ParticleGeometryResource;
class WignerDResources;

/**
 * @brief Maximum order and number of Gauss-Legendre quadrature points for a
 * converged T-matrix.
 */
class ConvergedSizeResources : public Resource {

  /*! @brief Grant access to @f$m=0@f$ computation task. */
  friend class TMatrixM0Task;

private:
  /*! @brief Iteration number. */
  uint_fast32_t _iteration;

  /*! @brief Maximum order, @f$n_{max}@f$. */
  uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  uint_fast32_t _ngauss;

  /*! @brief Gauss-Legendre quadrature points for the converged T-matrix. */
  const GaussBasedResources *_quadrature_points;

  /*! @brief Geometrical variables for the converged T-matrix. */
  const ParticleGeometryResource *_geometry;

  /*! @brief Wigner D functions for the converged T-matrix. */
  const WignerDResources *_wigner_d;

  /*! @brief Flag signaling whether or not the T-matrix caluclation is
   *  converged. */
  bool _is_converged;

  /*! @brief Relative difference for the scattering factor. */
  float_type _dQs;

  /*! @brief Relative difference for the extinction factor. */
  float_type _dQe;

public:
  /**
   * @brief Empty constructor.
   */
  inline ConvergedSizeResources()
      : _iteration(0), _nmax(0), _ngauss(0), _quadrature_points(nullptr),
        _geometry(nullptr), _wigner_d(nullptr), _is_converged(false), _dQs(-1.),
        _dQe(-1.) {}

  /**
   * @brief Get the size in memory of a hypothetical ConvergedSizeResources
   * object with the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(ConvergedSizeResources);
    return size;
  }

  /**
   * @brief Get the maximum order.
   *
   * @return Maximum order, @f$n_{max}@f$.
   */
  inline uint_fast32_t get_nmax() const { return _nmax; }

  /**
   * @brief Get the number of Gauss-Legendre quadrature points.
   *
   * @return Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$.
   */
  inline uint_fast32_t get_ngauss() const { return _ngauss; }

  /**
   * @brief Get the Gauss-Legendre quadrature points for the converged T-matrix.
   *
   * @return Gauss-Legendre quadrature points.
   */
  inline const GaussBasedResources *get_quadrature_points() const {
    return _quadrature_points;
  }

  /**
   * @brief Get the geometrical variables for the converged T-matrix.
   *
   * @return Geometrical variables.
   */
  inline const ParticleGeometryResource *get_geometry() const {
    return _geometry;
  }

  /**
   * @brief Get the Wigner D functions for the converged T-matrix.
   *
   * @return Wigner D functions.
   */
  inline const WignerDResources *get_wigner_d() const { return _wigner_d; }

  /**
   * @brief Is the T-matrix calculation converged?
   *
   * @return True if the calculation is converged.
   */
  inline bool is_converged() const { return _is_converged; }

  /**
   * @brief Get the relative difference for the scattering factor.
   *
   * @return Relative difference for the scattering factor.
   */
  inline float_type get_dQs() const { return _dQs; }

  /**
   * @brief Get the relative difference for the extinction factor.
   *
   * @return Relative difference for the extinction factor.
   */
  inline float_type get_dQe() const { return _dQe; }
};

#endif // CONVERGEDSIZERESOURCES_HPP
