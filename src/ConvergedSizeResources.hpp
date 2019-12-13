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
  /*! @brief Maximum order, @f$n_{max}@f$. */
  uint_fast32_t _nmax;

  /*! @brief Number of Gauss-Legendre quadrature points, @f$n_{GL}@f$. */
  uint_fast32_t _ngauss;

  /*! @brief Gauss-Legendre quadrature points for the converged T-matrix. */
  const GaussBasedResources *_quadrature_points;

  /*! @brief Geometrical variables for the converged T-matrix. */
  const ParticleGeometryResource *_geometry;

  /*! @brief Interaction variables for the converged T-matrix. */
  const InteractionResource *_interaction;

  /*! @brief Wigner D functions for the converged T-matrix. */
  const WignerDResources *_wigner_d;

  /*! @brief Flag signaling whether or not the T-matrix caluclation is
   *  converged. */
  bool _is_converged;

public:
  /**
   * @brief Empty constructor.
   */
  inline ConvergedSizeResources()
      : _nmax(0), _ngauss(0), _quadrature_points(nullptr), _geometry(nullptr),
        _interaction(nullptr), _wigner_d(nullptr), _is_converged(false) {}

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
   * @brief Get the interaction variables for the converged T-matrix.
   *
   * @return Interaction variables.
   */
  inline const InteractionResource *get_interaction() const {
    return _interaction;
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
};

#endif // CONVERGEDSIZERESOURCES_HPP
