/**
 * @file ShapeAveragingTask.hpp
 *
 * @brief Task that computes the shape distribution average of the absorption
 * coefficients.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SHAPEAVERAGINGTASK_HPP
#define SHAPEAVERAGINGTASK_HPP

#include "AbsorptionCoefficientTask.hpp"
#include "Configuration.hpp"
#include "QuickSchedWrapper.hpp"
#include "ShapeDistribution.hpp"

/**
 * @brief Task that computes the shape distribution average of the absorption
 * coefficients.
 */
class ShapeAveragingTask : public Task {
private:
  /*! @brief Shape distribution. */
  const ShapeDistribution &_shape_distribution;

  /*! @brief Input absorption coefficients. */
  std::vector<AbsorptionCoefficientResult *> _input_coefficients;

  /*! @brief Output (averaged) absorption coefficients. */
  AbsorptionCoefficientResult &_output_coefficients;

public:
  /**
   * @brief Constructor.
   *
   * @param shape_distribution Shape distribution.
   * @param output_coefficients Output coefficients.
   */
  inline ShapeAveragingTask(const ShapeDistribution &shape_distribution,
                            AbsorptionCoefficientResult &output_coefficients)
      : _shape_distribution(shape_distribution),
        _input_coefficients(shape_distribution.get_number_of_points(), nullptr),
        _output_coefficients(output_coefficients) {}

  virtual ~ShapeAveragingTask() {}

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _output_coefficients, true);
  }

  /**
   * @brief Add input coefficients for the calculation.
   *
   * @param quicksched QuickSched library.
   * @param ishape Index of the shape for which the coefficients are computed.
   * @param input_coefficient Input coefficients.
   */
  inline void
  add_input_coefficient(QuickSched &quicksched, const uint_fast32_t ishape,
                        AbsorptionCoefficientResult *input_coefficient) {
    ctm_assert(ishape < _input_coefficients.size());
    _input_coefficients[ishape] = input_coefficient;
    quicksched.link_task_and_resource(*this, *input_coefficient, false);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _output_coefficients._Qabs.size();
    const uint_fast32_t nshape = _shape_distribution.get_number_of_points();

    // make sure the average values are set to 0
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qabs[itheta] = 0.;
      _output_coefficients._Qabspol[itheta] = 0.;
    }

    // compute the nominator and denominator in the expression for the
    // average
    float_type norm = 0.;
    for (uint_fast32_t ishape = 0; ishape < nshape; ++ishape) {

      // some sanity checks
      ctm_assert(_input_coefficients[ishape] != nullptr);
      ctm_assert(_input_coefficients[ishape]->_Qabs.size() == ntheta);

      // get the value of the shape distribution at this evaluation point
      const float_type weight = _shape_distribution.get_weight(ishape);
      // add it to the norm (denominator in expression)
      norm += weight;
      // add the contributions to the absorption coefficients from this shape
      for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
        _output_coefficients._Qabs[itheta] +=
            weight * _input_coefficients[ishape]->_Qabs[itheta];
        _output_coefficients._Qabspol[itheta] +=
            weight * _input_coefficients[ishape]->_Qabspol[itheta];
      }
    }

    // normalise the average quantities
    const float_type norm_inv = 1. / norm;
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qabs[itheta] *= norm_inv;
      _output_coefficients._Qabspol[itheta] *= norm_inv;
    }
  }
};

#endif // SHAPEAVERAGINGTASK_HPP
