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

#include "AbsorptionCoefficientTask.hpp"
#include "AlignmentAverageTask.hpp"
#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "DraineDustProperties.hpp"
#include "DraineHensleyShapeDistribution.hpp"
#include "GaussBasedResources.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "NBasedResources.hpp"
#include "OrientationDistributionResource.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "ShapeAveragingTask.hpp"
#include "ShapeDistribution.hpp"
#include "SingleShapeShapeDistribution.hpp"
#include "TMatrixResource.hpp"
#include "Utilities.hpp"
#include "WignerDResources.hpp"

#include <cinttypes>
#include <cstddef>
#include <vector>

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

  /*! @brief Orientation distribution for non-aligned particles. */
  OrientationDistribution *_non_aligned_orientation_distribution;

  /*! @brief Orientation distribution for oblate aligned particles. */
  OrientationDistribution *_oblate_aligned_orientation_distribution;

  /*! @brief Orientation distribution for prolate aligned particles. */
  OrientationDistribution *_prolate_aligned_orientation_distribution;

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
   * @param aligned_orientation_distribution_type Type of orientation
   * distribution for aligned particles. Needs to be documented properly.
   */
  inline TaskManager(const uint_fast32_t minimum_order,
                     const uint_fast32_t maximum_order,
                     const uint_fast32_t gauss_legendre_factor,
                     const float_type tolerance,
                     const size_t maximum_memory_usage,
                     const int_fast32_t shape_distribution_type,
                     const int_fast32_t aligned_orientation_distribution_type)
      : _minimum_order(minimum_order), _maximum_order(maximum_order),
        _gauss_legendre_factor(gauss_legendre_factor), _tolerance(tolerance),
        _maximum_memory_usage(maximum_memory_usage),
        _shape_distribution(nullptr) {

    if (shape_distribution_type == 0) {
      _shape_distribution = new ShapeDistribution();
      _shape_distribution->evaluate(100u);
    } else if (shape_distribution_type == 1) {
      _shape_distribution = new DraineHensleyShapeDistribution(100u);
    } else if (shape_distribution_type == 2) {
      _shape_distribution = new SingleShapeShapeDistribution(1.00001);
    }

    _non_aligned_orientation_distribution =
        new OrientationDistribution(2 * _maximum_order);
    _non_aligned_orientation_distribution->initialise();
    if (aligned_orientation_distribution_type == 0) {
      _oblate_aligned_orientation_distribution =
          new DavisGreensteinOrientationDistribution(2 * _maximum_order, 2.);
      _prolate_aligned_orientation_distribution =
          new DavisGreensteinOrientationDistribution(2 * _maximum_order, 0.5);
    } else if (aligned_orientation_distribution_type == 1) {
      // note that we hardcode the cos2beta values for now:
      //  - oblate: 3/5
      //  - prolate: 1/5
      _oblate_aligned_orientation_distribution =
          new MishchenkoOrientationDistribution(2 * _maximum_order, 0.6);
      _prolate_aligned_orientation_distribution =
          new MishchenkoOrientationDistribution(2 * _maximum_order, 0.2);
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~TaskManager() {
    delete _shape_distribution;
    delete _non_aligned_orientation_distribution;
    delete _oblate_aligned_orientation_distribution;
    delete _prolate_aligned_orientation_distribution;
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
   * @param theta Grid of absorption coefficient angles (in radians).
   * @param ngauss Number of Gauss-Legendre quadrature points to use to
   * compute directional averages.
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
   */
  inline void generate_tasks(const std::vector<float_type> &theta,
                             const uint_fast32_t ngauss, QuickSched &quicksched,
                             std::vector<Task *> &tasks,
                             std::vector<Resource *> &resources,
                             std::vector<Result *> &results,
                             TMatrixAuxiliarySpaceManager *&space_manager) {

    // sort the input arrays
    std::sort(_sizes.begin(), _sizes.end());
    std::sort(_wavelengths.begin(), _wavelengths.end());

    // get the dimensions of the computational grid
    const uint_fast32_t number_of_compositions = _compositions.size();
    const uint_fast32_t number_of_sizes = _sizes.size();
    const uint_fast32_t number_of_wavelengths = _wavelengths.size();
    const uint_fast32_t number_of_shapes =
        _shape_distribution->get_number_of_points();
    // number of angles to use to compute absorption coefficient grid
    const uint_fast32_t number_of_angles = theta.size();

    // grand total of T-matrices we need to compute
    const uint_fast32_t total_number_of_interactions =
        number_of_compositions * number_of_sizes * number_of_wavelengths;
    const uint_fast32_t total_number_of_Tmatrices =
        total_number_of_interactions * number_of_shapes;
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
    nbased_resources->link_resources(quicksched);

    //  - orientation distribution resources
    OrientationDistributionResource *random_orientation =
        new OrientationDistributionResource(
            _non_aligned_orientation_distribution);
    OrientationDistributionResource *oblate_orientation =
        new OrientationDistributionResource(
            _oblate_aligned_orientation_distribution);
    OrientationDistributionResource *prolate_orientation =
        new OrientationDistributionResource(
            _prolate_aligned_orientation_distribution);
    quicksched.register_resource(*random_orientation);
    quicksched.register_resource(*oblate_orientation);
    quicksched.register_resource(*prolate_orientation);
    resources.resize(3, nullptr);
    resources[0] = random_orientation;
    resources[1] = oblate_orientation;
    resources[2] = prolate_orientation;

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
      this_quadrature_points->link_resources(quicksched);
      quadrature_points[i] = this_quadrature_points;
      tasks[quadrature_points_offset + i] = this_quadrature_points;
    }

    //  - absorption coefficient grid
    memory_used +=
        AbsorptionCoefficientGrid::get_memory_size(number_of_angles, ngauss);
    ctm_assert(memory_used <= _maximum_memory_usage);
    AbsorptionCoefficientGrid *grid =
        new AbsorptionCoefficientGrid(number_of_angles, &theta[0], ngauss);
    quicksched.register_resource(*grid);
    quicksched.register_task(*grid);
    grid->link_resources(quicksched);
    tasks.push_back(grid);

    //  - special Wigner D functions
    memory_used +=
        SpecialWignerDResources::get_memory_size(_maximum_order, *grid);
    ctm_assert(memory_used <= _maximum_memory_usage);
    SpecialWignerDResources *special_wigner =
        new SpecialWignerDResources(_maximum_order, *grid);
    quicksched.register_resource(*special_wigner);
    quicksched.register_task(*special_wigner);
    special_wigner->link_resources(quicksched);
    quicksched.link_tasks(*grid, *special_wigner);
    tasks.push_back(special_wigner);

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
        ParticleGeometryResource *this_geometry =
            new ParticleGeometryResource(_shape_distribution->get_shape(is),
                                         this_ngauss, *quadrature_points[ig]);
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

    // make sure we have enough memory left to store results
    memory_used -=
        total_number_of_Tmatrices *
        AbsorptionCoefficientResult::get_memory_size(number_of_angles);

    const DraineDustProperties dust_properties;
    // step 4: loop over all parameter values and set up parameter specific
    // tasks
    // first: figure out how much space is left for interaction and T-matrix
    // resources
    const uint_fast32_t memory_left = _maximum_memory_usage - memory_used;
    const uint_fast32_t tmatrix_memory_requirement =
        InteractionResource::get_memory_size(_maximum_order, maximum_ngauss) +
        2 * TMatrixResource::get_memory_size(_maximum_order);
    // compute the number of T-matrices that we can store simultaneously
    // if this number is larger than what we need, we only allocate what
    // we need
    const uint_fast32_t number_of_tmatrices = std::min(
        memory_left / tmatrix_memory_requirement, total_number_of_Tmatrices);
    // make sure the number we can allocate is larger than the number of
    // shapes (not strictly necessary, but we do this for now)
    ctm_assert(number_of_tmatrices >= number_of_shapes);

    // we are done using memory: output the total
    memory_used += number_of_tmatrices * tmatrix_memory_requirement;
    ctm_warning("Total memory usage: %s",
                Utilities::human_readable_bytes(memory_used).c_str());
    ctm_warning("Space to store %" PRIuFAST32 " T-matrices.",
                number_of_tmatrices);

    // now actually allocate the resources
    std::vector<InteractionResource *> interaction_resources(
        number_of_tmatrices, nullptr);
    std::vector<TMatrixResource *> tmatrices(2 * number_of_tmatrices, nullptr);
    const uint_fast32_t tmatrix_offset = resources.size();
    resources.resize(tmatrix_offset + 3 * number_of_tmatrices, nullptr);
    for (uint_fast32_t i = 0; i < number_of_tmatrices; ++i) {
      interaction_resources[i] =
          new InteractionResource(_maximum_order, maximum_ngauss);
      quicksched.register_resource(*interaction_resources[i]);
      resources[tmatrix_offset + 3 * i] = interaction_resources[i];

      tmatrices[2 * i] = new TMatrixResource(_maximum_order);
      resources[tmatrix_offset + 3 * i + 1] = tmatrices[2 * i];
      quicksched.register_resource(*tmatrices[2 * i]);
      for (uint_fast32_t m = 0; m < _maximum_order + 1; ++m) {
        quicksched.register_resource(tmatrices[2 * i]->get_m_resource(m),
                                     tmatrices[2 * i]);
      }
      tmatrices[2 * i + 1] = new TMatrixResource(_maximum_order);
      resources[tmatrix_offset + 3 * i + 2] = tmatrices[2 * i + 1];
      quicksched.register_resource(*tmatrices[2 * i + 1]);
      for (uint_fast32_t m = 0; m < _maximum_order + 1; ++m) {
        quicksched.register_resource(tmatrices[2 * i + 1]->get_m_resource(m),
                                     tmatrices[2 * i + 1]);
      }
    }

    // memory management:
    // we only have N = number_of_tmatrices available interaction resources
    // and T-matrix resources (2 of the latter, but that doesn't really matter
    // here). This means that we need to reuse these resources if we need more.
    // This is of course only allowed if we no longer need these resources, i.e.
    // if the AbsorptionCoefficientTask for that specific T-matrix has
    // finished.
    // We keep a running index of the last used T-matrix resources, and if that
    // overflows, we start reusing values. We also store all the tasks that
    // last used the resource and create a task dependency between the old
    // task that used the resource and the new task.
    std::vector<Task *> unlock_tasks(number_of_tmatrices, nullptr);
    uint_fast32_t running_resource_index = 0;

    // big loop over all parameter values
    // for each parameter value we allocate interaction variables and a
    // control resource
    std::vector<InteractionVariables *> interaction_variables(
        total_number_of_interactions, nullptr);
    uint_fast32_t interaction_index_new = 0;
    std::vector<ConvergedSizeResources *> converged_size_resources(
        total_number_of_Tmatrices, nullptr);
    uint_fast32_t converged_size_index = 0;
    std::vector<AbsorptionCoefficientResult *> unaveraged_results(
        total_number_of_Tmatrices, nullptr);
    uint_fast32_t unaveraged_result_index = 0;
    const uint_fast32_t resource_offset = resources.size();
    resources.resize(resource_offset + total_number_of_interactions +
                         2 * total_number_of_Tmatrices,
                     nullptr);
    uint_fast32_t resource_index = resource_offset;
    results.resize(total_number_of_interactions, nullptr);
    uint_fast32_t result_index = 0;
    const uint_fast32_t task_offset = tasks.size();
    // for each T-matrix,
    // we need to add interaction and m=0 tasks for each quadrature point
    // we need to add m=/=0 tasks for each order
    const uint_fast32_t tasks_per_Tmatrix =
        2 * number_of_quadrature_tasks + _maximum_order + 5;
    tasks.resize(task_offset + total_number_of_interactions +
                     total_number_of_Tmatrices * tasks_per_Tmatrix,
                 nullptr);
    uint_fast32_t task_index = task_offset;
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

          InteractionVariables *this_interaction_variables =
              new InteractionVariables(particle_size, wavelength,
                                       refractive_index);
          quicksched.register_resource(*this_interaction_variables);
          interaction_variables[interaction_index_new] =
              this_interaction_variables;
          ++interaction_index_new;
          resources[resource_index] = this_interaction_variables;
          ++resource_index;

          AbsorptionCoefficientResult *this_result =
              new AbsorptionCoefficientResult(grain_type, particle_size,
                                              wavelength, number_of_angles);
          quicksched.register_resource(*this_result);
          results[result_index] = this_result;
          ++result_index;

          ShapeAveragingTask *averaging_task =
              new ShapeAveragingTask(*_shape_distribution, *this_result);
          quicksched.register_task(*averaging_task);
          averaging_task->link_resources(quicksched);
          tasks[task_index] = averaging_task;
          ++task_index;

          // loop over all shapes
          for (uint_fast32_t ishape = 0; ishape < number_of_shapes; ++ishape) {

            AbsorptionCoefficientResult *this_unaveraged_result =
                new AbsorptionCoefficientResult(grain_type, particle_size,
                                                wavelength, number_of_angles);
            quicksched.register_resource(*this_unaveraged_result);
            unaveraged_results[unaveraged_result_index] =
                this_unaveraged_result;
            ++unaveraged_result_index;
            resources[resource_index] = this_unaveraged_result;
            ++resource_index;

            OrientationDistributionResource *this_orientation;
            if (particle_size > 1.e-5) {
              if (_shape_distribution->get_shape(ishape) > 1.) {
                this_orientation = oblate_orientation;
              } else {
                this_orientation = prolate_orientation;
              }
            } else {
              this_orientation = random_orientation;
            }

            // allocate the corresponding interaction variables and control
            // object
            ConvergedSizeResources *this_converged_size =
                new ConvergedSizeResources();
            quicksched.register_resource(*this_converged_size);
            converged_size_resources[converged_size_index] =
                this_converged_size;
            ++converged_size_index;
            resources[resource_index] = this_converged_size;
            ++resource_index;

            InteractionResource *this_interaction_resource =
                interaction_resources[running_resource_index];
            TMatrixResource *this_single_Tmatrix =
                tmatrices[2 * running_resource_index];
            TMatrixResource *this_ensemble_Tmatrix =
                tmatrices[2 * running_resource_index + 1];

            // we need to store the previous m=0 task to set up a dependency
            TMatrixM0Task *previous_m0task = nullptr;
            // loop over all orders
            for (uint_fast32_t ig = 0; ig < number_of_quadrature_tasks; ++ig) {
              const uint_fast32_t this_order = _minimum_order + ig;
              const uint_fast32_t this_ngauss =
                  minimum_ngauss + ig * _gauss_legendre_factor;
              const GaussBasedResources &this_quadrature_points =
                  *quadrature_points[ig];
              const WignerDResources &this_wigner = *wignerdm0[ig];
              const ParticleGeometryResource &this_geometry =
                  *geometries[ishape * number_of_quadrature_tasks + ig];
              // set up the interaction task
              InteractionTask *this_interaction = new InteractionTask(
                  this_order, this_ngauss, this_geometry, *this_converged_size,
                  *this_interaction_variables, *this_interaction_resource);
              quicksched.register_task(*this_interaction);
              this_interaction->link_resources(quicksched);
              quicksched.link_tasks(this_geometry, *this_interaction);
              if (unlock_tasks[running_resource_index] != nullptr) {
                quicksched.link_tasks(*unlock_tasks[running_resource_index],
                                      *this_interaction);
              }
              tasks[task_index] = this_interaction;
              ++task_index;

              // set up the m=0 task
              TMatrixM0Task *this_m0task = new TMatrixM0Task(
                  _tolerance, this_order, this_ngauss, *nbased_resources,
                  this_quadrature_points, this_geometry,
                  *this_interaction_variables, *this_interaction_resource,
                  this_wigner, *space_manager, *this_single_Tmatrix,
                  *this_converged_size, this_single_Tmatrix->get_m_resource(0));
              quicksched.register_task(*this_m0task);
              this_m0task->link_resources(quicksched);
              quicksched.link_tasks(*nbased_resources, *this_m0task);
              quicksched.link_tasks(this_wigner, *this_m0task);
              quicksched.link_tasks(*this_interaction, *this_m0task);
              // link the previous task as dependency, if it exists
              if (previous_m0task != nullptr) {
                quicksched.link_tasks(*previous_m0task, *this_interaction);
              } else {
                if (unlock_tasks[running_resource_index] != nullptr) {
                  quicksched.link_tasks(*unlock_tasks[running_resource_index],
                                        *this_m0task);
                }
              }
              tasks[task_index] = this_m0task;
              ++task_index;
              previous_m0task = this_m0task;
            }

            AlignmentAverageTask *alignment_task = new AlignmentAverageTask(
                *this_orientation, *this_single_Tmatrix,
                *this_ensemble_Tmatrix);
            quicksched.register_task(*alignment_task);
            alignment_task->link_resources(quicksched);
            tasks[task_index] = alignment_task;
            ++task_index;

            AbsorptionCoefficientTask *absorption_task =
                new AbsorptionCoefficientTask(
                    *grid, *this_interaction_variables, *this_ensemble_Tmatrix,
                    *nbased_resources, *special_wigner,
                    *this_unaveraged_result);
            quicksched.register_task(*absorption_task);
            absorption_task->link_resources(quicksched);
            quicksched.link_tasks(*alignment_task, *absorption_task);
            quicksched.link_tasks(*special_wigner, *absorption_task);
            tasks[task_index] = absorption_task;
            ++task_index;

            averaging_task->add_input_coefficient(quicksched, ishape,
                                                  this_unaveraged_result);
            quicksched.link_tasks(*absorption_task, *averaging_task);

            ResetTMatrixResourceTask *reset_task1 =
                new ResetTMatrixResourceTask(*this_single_Tmatrix);
            quicksched.register_task(*reset_task1);
            reset_task1->link_resources(quicksched);
            tasks[task_index] = reset_task1;
            ++task_index;
            quicksched.link_tasks(*alignment_task, *reset_task1);

            ResetTMatrixResourceTask *reset_task2 =
                new ResetTMatrixResourceTask(*this_ensemble_Tmatrix);
            quicksched.register_task(*reset_task2);
            reset_task2->link_resources(quicksched);
            tasks[task_index] = reset_task2;
            ++task_index;
            // this one is not stricly necessary but makes it easier to keep
            // track of unlock tasks
            quicksched.link_tasks(*reset_task1, *reset_task2);
            quicksched.link_tasks(*absorption_task, *reset_task2);

            ResetInteractionResourceTask *reset_task3 =
                new ResetInteractionResourceTask(*this_interaction_resource);
            quicksched.register_task(*reset_task3);
            reset_task3->link_resources(quicksched);
            tasks[task_index] = reset_task3;
            ++task_index;
            // this one is not stricly necessary but makes it easier to keep
            // track of unlock tasks
            quicksched.link_tasks(*reset_task2, *reset_task3);

            // loop over all m values to set up m=/=0 tasks
            for (uint_fast32_t i = 0; i < _maximum_order; ++i) {
              TMatrixMAllTask *this_malltask = new TMatrixMAllTask(
                  i + 1, *nbased_resources, *this_converged_size,
                  *this_interaction_variables, *this_interaction_resource,
                  *space_manager, *this_single_Tmatrix,
                  this_single_Tmatrix->get_m_resource(i + 1));
              quicksched.register_task(*this_malltask);
              this_malltask->link_resources(quicksched);
              // link to the last m=0 task
              quicksched.link_tasks(*previous_m0task, *this_malltask);
              quicksched.link_tasks(*this_malltask, *alignment_task);
              tasks[task_index] = this_malltask;
              ++task_index;
            }

            unlock_tasks[running_resource_index] = reset_task3;
            ++running_resource_index;
            running_resource_index %= number_of_tmatrices;
          }
        }
      }
    }
  }
};

#endif // TASKMANAGER_HPP
