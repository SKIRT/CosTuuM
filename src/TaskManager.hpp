/**
 * @file TaskManager.hpp
 *
 * @brief Class that generates tasks and dependencies for a set of T-matrix
 * calculations.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TASKMANAGER_HPP
#define TASKMANAGER_HPP

#include "Configuration.hpp"
#include "DraineDustProperties.hpp"
#include "DraineHensleyShapeDistribution.hpp"
#include "GaussBasedResources.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "ShapeDistribution.hpp"
#include "TMatrixResource.hpp"
#include "Utilities.hpp"
#include "WignerDResources.hpp"

#include <cinttypes>
#include <cstddef>
#include <vector>

/**
 * @brief Interface for task-based results.
 */
class Result : public Resource {
private:
  /*! @brief Composition parameter value for the result. */
  const int_fast32_t _composition;

  /*! @brief Particle size parameter value for the result. */
  const float_type _size;

  /*! @brief Wavelength value for the result. */
  const float_type _wavelength;

  /*! @brief Type of result. */
  const int_fast32_t _type;

public:
  /**
   * @brief Constructor.
   *
   * @param composition Composition parameter value for the result.
   * @param size Particle size parameter value for the result (in m).
   * @param wavelength Wavelength value for the result (in m).
   * @param type Type of result.
   */
  inline Result(const int_fast32_t composition, const float_type size,
                const float_type wavelength, const int_fast32_t type)
      : _composition(composition), _size(size), _wavelength(wavelength),
        _type(type) {}
};

/**
 * @brief Class that generates tasks and dependencies for a set of T-matrix
 * calculations.
 */
class TaskManager {
private:
  /*! @brief Minimum order of spherical basis function expansions; this is the
   *  starting point for all iterative T-matrix calculations. */
  const uint_fast32_t _minimum_order;

  /*! @brief Maximum allowed order of spherical basis function expansions. */
  const uint_fast32_t _maximum_order;

  /*! @brief Gauss-Legendre factor. The minimum and maximum allowed number of
   *  Gauss-Legendre quadrature points corresponds to respectively the minimum
   *  and maximum order multiplied with this factor. */
  const uint_fast32_t _gauss_legendre_factor;

  /*! @brief Tolerance for T-matrix calculations. */
  const float_type _tolerance;

  /*! @brief Maximum allowed memory usage (in bytes) for memory intensive
   *  resources used by the calculation. This will be used to limit the number
   *  of T-matrices that can be stored and computed simultaneously. */
  const size_t _maximum_memory_usage;

  /*! @brief Shape distribution for the particle shapes. */
  ShapeDistribution *_shape_distribution;

  /*! @brief Requested particle compositions. */
  std::vector<int_fast32_t> _compositions;

  /*! @brief Requested particle sizes (in m). */
  std::vector<float_type> _sizes;

  /*! @brief Requested wavelengths (in m). */
  std::vector<float_type> _wavelengths;

public:
  /**
   * @brief Constructor.
   *
   * @param minimum_order Minimum order of spherical basis function expansions.
   * @param maximum_order Maximum order of spherical basis function expansions.
   * @param gauss_legendre_factor Multiplicative factor setting the number of
   * Gauss-Legendre quadrature points for an expansion of a given order when
   * multiplying it with that order.
   * @param tolerance Tolerance used to decide when a T-matrix calculation is
   * converged.
   * @param maximum_memory_usage Maximum allowed memory usage for the algorithm
   * (in bytes). This value is only taken into account when allocating task-
   * based structures and might be exceeded a bit during some steps of the
   * algorithm.
   * @param shape_distribution_type Type of shape distribution. Needs to be
   * documented properly.
   */
  inline TaskManager(const uint_fast32_t minimum_order,
                     const uint_fast32_t maximum_order,
                     const uint_fast32_t gauss_legendre_factor,
                     const float_type tolerance,
                     const size_t maximum_memory_usage,
                     int_fast32_t shape_distribution_type)
      : _minimum_order(minimum_order), _maximum_order(maximum_order),
        _gauss_legendre_factor(gauss_legendre_factor), _tolerance(tolerance),
        _maximum_memory_usage(maximum_memory_usage),
        _shape_distribution(nullptr) {

    if (shape_distribution_type == 0) {
      _shape_distribution = new ShapeDistribution();
    } else if (shape_distribution_type == 1) {
      _shape_distribution = new DraineHensleyShapeDistribution();
    }
  }

  /**
   * @brief Add the given particle composition to the list of compositions that
   * need to be computed.
   *
   * @param composition Composition type.
   */
  inline void add_composition(const int_fast32_t composition) {
    _compositions.push_back(composition);
  }

  /**
   * @brief Add the given particle size to the list of sizes that need to be
   * computed.
   *
   * @param size Particle size (equal volume sphere radius, in m).
   */
  inline void add_size(const float_type size) { _sizes.push_back(size); }

  /**
   * @brief Add the given wavelength to the list of wavelengths that need to be
   * computed.
   *
   * @param wavelength Wavelength (in m).
   */
  inline void add_wavelength(const float_type wavelength) {
    _wavelengths.push_back(wavelength);
  }

  /**
   * @brief Generate and link all the tasks required to compute the configured
   * grid of parameters.
   *
   * @param quicksched QuickSched library wrapper.
   * @param tasks List of tasks. Tasks need to be deleted by caller after the
   * computation finishes. This list also contains resources that act as a task.
   * @param resources List of resources. Resources need to be deleted by caller
   * after the computation finishes.
   * @param results Results of the computation. Note that these are also
   * resources. Results need to be deleted by caller after the computation
   * finishes (and after the results have been retrieved).
   * @param space_manager TMatrixAuxiliarySpaceManager that manages per thread
   * additional space for the T-matrix calculation. Needs to be deleted by
   * caller after the computation finishes.
   * @param tmatrices T-matrix resources. Need to be deleted by caller.
   * @param interaction_variables InteractionVariables. Need to be deleted by
   * caller.
   */
  inline void generate_tasks(
      QuickSched &quicksched, std::vector<Task *> &tasks,
      std::vector<Resource *> &resources, std::vector<Result *> &results,
      TMatrixAuxiliarySpaceManager *space_manager,
      std::vector<TMatrixResource *> &tmatrices,
      std::vector<InteractionVariables *> &interaction_variables) const {

    // get the dimensions of the computational grid
    const uint_fast32_t number_of_compositions = _compositions.size();
    const uint_fast32_t number_of_sizes = _sizes.size();
    const uint_fast32_t number_of_wavelengths = _wavelengths.size();
    // I guess the number of shape quadrature points should somehow be
    // determined from the distribution or be a parameter
    // But for now, we simply use a hard coded value...
    const uint_fast32_t number_of_shapes = 100u;

    // grand total of T-matrices we need to compute
    const uint_fast32_t total_number_of_Tmatrices =
        number_of_compositions * number_of_sizes * number_of_wavelengths *
        number_of_shapes;
    ctm_warning("Number of T-matrices that needs to be computed: %" PRIuFAST32,
                total_number_of_Tmatrices);

    // we need to keep track of memory to make sure we respect the user
    // defined limit
    uint_fast32_t memory_used = 0;

    // allocate auxiliary space per thread
    memory_used += TMatrixAuxiliarySpaceManager::get_memory_size(
        quicksched.get_number_of_threads(), _maximum_order);
    ctm_assert(memory_used <= _maximum_memory_usage);
    space_manager = new TMatrixAuxiliarySpaceManager(
        quicksched.get_number_of_threads(), _maximum_order);

    // step 1: set up the tasks that need to be done regardless of what the
    // parameter values are:

    //  - N based resources
    memory_used += NBasedResources::get_memory_size(_maximum_order);
    ctm_assert(memory_used <= _maximum_memory_usage);
    NBasedResources *nbased_resources = new NBasedResources(_maximum_order);
    tasks.push_back(nbased_resources);
    quicksched.register_resource(*nbased_resources);
    quicksched.register_task(*nbased_resources);

    //  - quadrature points
    const uint_fast32_t minimum_ngauss =
        _minimum_order * _gauss_legendre_factor;
    const uint_fast32_t maximum_ngauss =
        _maximum_order * _gauss_legendre_factor;
    const uint_fast32_t number_of_quadrature_tasks =
        _maximum_order - _minimum_order;
    std::vector<GaussBasedResources *> quadrature_points(
        number_of_quadrature_tasks, nullptr);
    const uint_fast32_t quadrature_points_offset = tasks.size();
    tasks.resize(quadrature_points_offset + number_of_quadrature_tasks,
                 nullptr);
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      memory_used += GaussBasedResources::get_memory_size(this_ngauss);
      ctm_assert(memory_used <= _maximum_memory_usage);
      GaussBasedResources *this_quadrature_points =
          new GaussBasedResources(this_ngauss);
      quicksched.register_resource(*this_quadrature_points);
      quicksched.register_task(*this_quadrature_points);
      quadrature_points[i] = this_quadrature_points;
      tasks[quadrature_points_offset + i] = this_quadrature_points;
    }

    //  - Wigner D functions
    std::vector<WignerDResources *> wignerdm0(number_of_quadrature_tasks,
                                              nullptr);
    const uint_fast32_t wignerdm0_offset = tasks.size();
    tasks.resize(wignerdm0_offset + number_of_quadrature_tasks, nullptr);
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_order = _minimum_order + i;
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      memory_used += WignerDResources::get_memory_size(this_order, this_ngauss);
      ctm_assert(memory_used <= _maximum_memory_usage);
      WignerDResources *this_wignerdm0 =
          new WignerDResources(this_order, this_ngauss, *quadrature_points[i]);
      quicksched.register_resource(*this_wignerdm0);
      quicksched.register_task(*this_wignerdm0);
      this_wignerdm0->link_resources(quicksched);
      quicksched.link_tasks(*quadrature_points[i], *this_wignerdm0);
      wignerdm0[i] = this_wignerdm0;
      tasks[wignerdm0_offset + i] = this_wignerdm0;
    }

    // step 2: compute shape quadrature points and generate shape based
    // resources that are shared between all parameter values
    std::vector<float_type> shape_values(number_of_shapes),
        shape_weights(number_of_shapes);
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab(
        number_of_shapes, _shape_distribution->get_minimum_axis_ratio(),
        _shape_distribution->get_maximum_axis_ratio(), shape_values,
        shape_weights);
    const uint_fast32_t number_of_geometry_tasks =
        number_of_shapes * number_of_quadrature_tasks;
    std::vector<ParticleGeometryResource *> geometries(number_of_geometry_tasks,
                                                       nullptr);
    const uint_fast32_t shape_offset = tasks.size();
    tasks.resize(shape_offset + number_of_geometry_tasks, nullptr);
    for (uint_fast32_t is = 0; is < number_of_shapes; ++is) {
      for (uint_fast32_t ig = 0; ig < number_of_quadrature_tasks; ++ig) {
        const uint_fast32_t this_ngauss =
            minimum_ngauss + ig * _gauss_legendre_factor;
        memory_used += ParticleGeometryResource::get_memory_size(this_ngauss);
        ctm_assert(memory_used <= _maximum_memory_usage);
        ParticleGeometryResource *this_geometry = new ParticleGeometryResource(
            shape_values[is], this_ngauss, *quadrature_points[ig]);
        quicksched.register_resource(*this_geometry);
        quicksched.register_task(*this_geometry);
        this_geometry->link_resources(quicksched);
        quicksched.link_tasks(*quadrature_points[ig], *this_geometry);
        // note that we store the particle geometries as follows:
        // rows (first index) correspond to different shapes
        // columns (second index) correspond to different quadrature points
        const uint_fast32_t index = is * number_of_quadrature_tasks + ig;
        geometries[index] = this_geometry;
        tasks[shape_offset + index] = this_geometry;
      }
    }

    const DraineDustProperties dust_properties;
    // step 4: loop over all parameter values and set up parameter specific
    // tasks
    // first: figure out how much space is left for interaction and T-matrix
    // resources
    const uint_fast32_t memory_left = _maximum_memory_usage - memory_used;
    const uint_fast32_t tmatrix_memory_requirement =
        InteractionResource::get_memory_size(_maximum_order, maximum_ngauss) +
        TMatrixResource::get_memory_size(_maximum_order);
    // compute the number of T-matrices that we can store simultaneously
    // if this number is larger than what we need, we only allocate what
    // we need
    const uint_fast32_t number_of_tmatrices = std::min(
        memory_left / tmatrix_memory_requirement, total_number_of_Tmatrices);
    // make sure the number we can allocate is larger than the number of
    // shapes (not strictly necessary, but we do this for now)
    ctm_assert(number_of_tmatrices >= number_of_shapes);

    // for now: enforce enough space to store all T-matrices
    ctm_assert(number_of_tmatrices == total_number_of_Tmatrices);

    // we are done using memory: output the total
    memory_used += number_of_tmatrices * tmatrix_memory_requirement;
    ctm_warning("Total memory usage: %s",
                Utilities::human_readable_bytes(memory_used).c_str());

    // now actually allocate the resources
    std::vector<InteractionResource *> interaction_resources(
        number_of_tmatrices, nullptr);
    tmatrices.resize(number_of_tmatrices, nullptr);
    const uint_fast32_t tmatrix_offset = resources.size();
    resources.resize(tmatrix_offset + number_of_tmatrices, nullptr);
    for (uint_fast32_t i = 0; i < number_of_tmatrices; ++i) {
      interaction_resources[i] =
          new InteractionResource(_maximum_order, maximum_ngauss);
      quicksched.register_resource(*interaction_resources[i]);
      resources[tmatrix_offset + i] = interaction_resources[i];

      tmatrices[i] = new TMatrixResource(_maximum_order);
      for (uint_fast32_t m = 0; m < _maximum_order + 1; ++m) {
        quicksched.register_resource(tmatrices[i]->get_m_resource(m));
      }
    }

    // big loop over all parameter values
    // for each parameter value we allocate interaction variables and a
    // control resource
    interaction_variables.resize(total_number_of_Tmatrices, nullptr);
    std::vector<ConvergedSizeResources *> converged_size_resources(
        total_number_of_Tmatrices, nullptr);
    const uint_fast32_t resource_offset = resources.size();
    resources.resize(resource_offset + total_number_of_Tmatrices, nullptr);
    const uint_fast32_t task_offset = tasks.size();
    // for each T-matrix,
    // we need to add interaction and m=0 tasks for each quadrature point
    // we need to add m=/=0 tasks for each order
    tasks.resize(task_offset +
                     total_number_of_Tmatrices *
                         (2 * number_of_quadrature_tasks + _maximum_order),
                 nullptr);
    // loop over compositions
    for (uint_fast32_t icomp = 0; icomp < number_of_compositions; ++icomp) {
      // get the grain type that corresponds to this composition
      const int_fast32_t grain_type = _compositions[icomp];
      // make sure the grain type exists
      ctm_assert(grain_type >= 0 && grain_type < NUMBER_OF_DUSTGRAINTYPES);
      // loop over all sizes
      for (uint_fast32_t isize = 0; isize < number_of_sizes; ++isize) {
        // cache the size
        const float_type particle_size = _sizes[isize];
        // loop over all wavelengths
        for (uint_fast32_t ilambda = 0; ilambda < number_of_wavelengths;
             ++ilambda) {
          // cache the wavelength
          const float_type wavelength = _wavelengths[ilambda];
          // get the refractive index for this grain type, particle size and
          // interaction wavelength
          const std::complex<float_type> refractive_index =
              dust_properties.get_refractive_index(wavelength, particle_size,
                                                   grain_type);
          // loop over all shapes
          for (uint_fast32_t ishape = 0; ishape < number_of_shapes; ++ishape) {
            // compute the index of the corresponding T-matrix
            const uint_fast32_t index =
                icomp * number_of_sizes * number_of_wavelengths *
                    number_of_shapes +
                isize * number_of_wavelengths * number_of_shapes +
                ilambda * number_of_shapes + ishape;

            // allocate the corresponding interaction variables and control
            // object
            interaction_variables[index] = new InteractionVariables(
                particle_size, wavelength, refractive_index);
            converged_size_resources[index] = new ConvergedSizeResources();
            quicksched.register_resource(*converged_size_resources[index]);
            resources[resource_offset + index] =
                converged_size_resources[index];
            // we need to store the previous m=0 task to set up a dependency
            TMatrixM0Task *previous_m0task = nullptr;
            // loop over all orders
            for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
              const uint_fast32_t this_order = _minimum_order + i;
              const uint_fast32_t this_ngauss =
                  minimum_ngauss + i * _gauss_legendre_factor;
              const GaussBasedResources &this_quadrature_points =
                  *quadrature_points[i];
              const WignerDResources &this_wigner = *wignerdm0[i];
              const ParticleGeometryResource &this_geometry =
                  *geometries[ishape * number_of_quadrature_tasks + i];
              // set up the interaction task
              InteractionTask *this_interaction = new InteractionTask(
                  this_order, this_ngauss, this_geometry,
                  *converged_size_resources[index],
                  *interaction_variables[index], *interaction_resources[index]);
              quicksched.register_task(*this_interaction);
              this_interaction->link_resources(quicksched);
              quicksched.link_tasks(this_geometry, *this_interaction);
              tasks[task_offset +
                    index * (2 * number_of_quadrature_tasks + _maximum_order) +
                    2 * i] = this_interaction;

              // set up the m=0 task
              TMatrixM0Task *this_m0task = new TMatrixM0Task(
                  _tolerance, this_order, this_ngauss, *nbased_resources,
                  this_quadrature_points, this_geometry,
                  *interaction_variables[index], *interaction_resources[index],
                  this_wigner, *space_manager, *tmatrices[index],
                  *converged_size_resources[index],
                  tmatrices[index]->get_m_resource(0));
              quicksched.register_task(*this_m0task);
              this_m0task->link_resources(quicksched);
              quicksched.link_tasks(*nbased_resources, *this_m0task);
              quicksched.link_tasks(this_wigner, *this_m0task);
              quicksched.link_tasks(*this_interaction, *this_m0task);
              // link the previous task as dependency, if it exists
              if (previous_m0task != nullptr) {
                quicksched.link_tasks(*previous_m0task, *this_interaction);
                quicksched.link_tasks(*previous_m0task, *this_m0task);
              }
              tasks[task_offset +
                    index * (2 * number_of_quadrature_tasks + _maximum_order) +
                    2 * i + 1] = this_m0task;
              previous_m0task = this_m0task;
            }

            // loop over all m values to set up m=/=0 tasks
            for (uint_fast32_t i = 0; i < _maximum_order; ++i) {
              TMatrixMAllTask *this_malltask = new TMatrixMAllTask(
                  i + 1, *nbased_resources, *converged_size_resources[index],
                  *interaction_variables[index], *space_manager,
                  *tmatrices[index], tmatrices[index]->get_m_resource(i + 1));
              quicksched.register_task(*this_malltask);
              this_malltask->link_resources(quicksched);
              // link to the last m=0 task
              quicksched.link_tasks(*previous_m0task, *this_malltask);
              tasks[task_offset +
                    index * (2 * number_of_quadrature_tasks + _maximum_order) +
                    2 * number_of_quadrature_tasks + i] = this_malltask;
            }
          }
        }
      }
    }
  }
};

#endif // TASKMANAGER_HPP
