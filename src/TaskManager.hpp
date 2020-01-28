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
#include "AlignmentDistribution.hpp"
#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "DraineHensleyShapeDistribution.hpp"
#include "DustProperties.hpp"
#include "GaussBasedResources.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "NBasedResources.hpp"
#include "ParticleGeometryResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "ResultKey.hpp"
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
  const ShapeDistribution &_shape_distribution;

  /*! @brief Alignment distribution for the dust grains. */
  const AlignmentDistribution &_alignment_distribution;

  /*! @brief Dust grain properties. */
  const DustProperties &_dust_properties;

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
   * @param shape_distribution Shape distribution that specifies how the shapes
   * of dust grains are distributed.
   * @param alignment_distribution Alignment distribution that specifies how
   * dust grains align with the magnetic field.
   * @param dust_properties Dust grain properties.
   */
  inline TaskManager(const uint_fast32_t minimum_order,
                     const uint_fast32_t maximum_order,
                     const uint_fast32_t gauss_legendre_factor,
                     const float_type tolerance,
                     const size_t maximum_memory_usage,
                     const ShapeDistribution &shape_distribution,
                     const AlignmentDistribution &alignment_distribution,
                     const DustProperties &dust_properties)
      : _minimum_order(minimum_order), _maximum_order(maximum_order),
        _gauss_legendre_factor(gauss_legendre_factor), _tolerance(tolerance),
        _maximum_memory_usage(maximum_memory_usage),
        _shape_distribution(shape_distribution),
        _alignment_distribution(alignment_distribution),
        _dust_properties(dust_properties) {}

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
   * @brief Account for the given additional memory usage.
   *
   * @param size Additional memory usage size (in bytes).
   * @param label Label to mention in log file and error messages.
   * @param memory_log Memory log file to write to (or nullptr for no log).
   * @param cumulative_size Cumulative size of all allocations.
   */
  inline void add_memory_allocation(const size_t size, const std::string label,
                                    std::ofstream *memory_log,
                                    size_t &cumulative_size) const {

    if (memory_log) {
      *memory_log << label << "\t" << size << "\n";
      memory_log->flush();
    }
    cumulative_size += size;
    ctm_assert_message(cumulative_size < _maximum_memory_usage, "%s, %zu bytes",
                       label.c_str(), size);
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
   * @param result_key Object that helps obtain the index for a certain
   * parameter result from the results array when looping over the input
   * parameters in an ordered fashion. Needs to be deleted by caller after the
   * results have been retrieved.
   * @param results Results of the computation. Note that these are also
   * resources. Results need to be deleted by caller after the computation
   * finishes (and after the results have been retrieved).
   * @param space_manager TMatrixAuxiliarySpaceManager that manages per thread
   * additional space for the T-matrix calculation. Needs to be deleted by
   * caller after the computation finishes.
   * @param do_extinction Compute extinction coefficients?
   * @param do_absorption Compute absorption coefficients?
   * @param do_scattering Compute scattering matrices?
   * @param account_for_scattering Subtract the directionally averaged
   * scattering cross sections from the absorption coefficients (unstable)?
   * @param verbose Output diagnostic warnings?
   * @param write_memory_log Write a log file with the memory allocations?
   * @param memory_log_file_name Name of the memory log file.
   */
  inline void
  generate_tasks(const std::vector<float_type> &theta,
                 const uint_fast32_t ngauss, QuickSched &quicksched,
                 std::vector<Task *> &tasks, std::vector<Resource *> &resources,
                 ResultKey *&result_key, std::vector<Result *> &results,
                 TMatrixAuxiliarySpaceManager *&space_manager,
                 const bool do_extinction, const bool do_absorption,
                 const bool do_scattering, const bool account_for_scattering,
                 const bool verbose = false,
                 const bool write_memory_log = false,
                 const std::string memory_log_file_name = "") {

    std::ofstream *memory_log_file = nullptr;
    if (write_memory_log) {
      memory_log_file = new std::ofstream(memory_log_file_name);
      *memory_log_file << "# label\tallocation (bytes)\n";
      memory_log_file->flush();
    }

    // sort the input arrays
    std::sort(_sizes.begin(), _sizes.end());
    std::sort(_wavelengths.begin(), _wavelengths.end());

    uint_fast32_t number_of_results = 0;
    if (do_extinction) {
      ++number_of_results;
    }
    if (do_absorption) {
      ++number_of_results;
    }
    if (do_scattering) {
      ++number_of_results;
    }
    result_key =
        new ResultKey(_compositions, _sizes, _wavelengths, number_of_results);

    // get the dimensions of the computational grid
    const uint_fast32_t number_of_compositions = _compositions.size();
    const uint_fast32_t number_of_sizes = _sizes.size();
    const uint_fast32_t number_of_wavelengths = _wavelengths.size();
    const uint_fast32_t number_of_shapes =
        _shape_distribution.get_number_of_points();
    // number of angles to use to compute absorption coefficient grid
    const uint_fast32_t number_of_angles = theta.size();

    // grand total of T-matrices we need to compute
    const uint_fast32_t total_number_of_interactions =
        number_of_compositions * number_of_sizes * number_of_wavelengths;
    const uint_fast32_t total_number_of_Tmatrices =
        total_number_of_interactions * number_of_shapes;
    if (verbose) {
      ctm_warning(
          "Number of T-matrices that needs to be computed: %" PRIuFAST32,
          total_number_of_Tmatrices);
    }

    // we need to keep track of memory to make sure we respect the user
    // defined limit
    size_t memory_used = 0;

    // allocate auxiliary space per thread
    add_memory_allocation(
        TMatrixAuxiliarySpaceManager::get_memory_size(
            quicksched.get_number_of_threads(), _maximum_order),
        "TMatrixAuxiliarySpaceManager", memory_log_file, memory_used);
    space_manager = new TMatrixAuxiliarySpaceManager(
        quicksched.get_number_of_threads(), _maximum_order);

    // step 1: set up the tasks that need to be done regardless of what the
    // parameter values are:

    //  - N based resources
    add_memory_allocation(NBasedResources::get_memory_size(_maximum_order),
                          "NBasedResources", memory_log_file, memory_used);
    NBasedResources *nbased_resources = new NBasedResources(_maximum_order);
    tasks.push_back(nbased_resources);
    quicksched.register_resource(*nbased_resources);
    quicksched.register_task(*nbased_resources);
    nbased_resources->link_resources(quicksched);

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
      add_memory_allocation(GaussBasedResources::get_memory_size(this_ngauss),
                            "GaussBasedResources", memory_log_file,
                            memory_used);
      GaussBasedResources *this_quadrature_points =
          new GaussBasedResources(this_ngauss);
      quicksched.register_resource(*this_quadrature_points);
      quicksched.register_task(*this_quadrature_points);
      this_quadrature_points->link_resources(quicksched);
      quadrature_points[i] = this_quadrature_points;
      tasks[quadrature_points_offset + i] = this_quadrature_points;
    }

    //  - extinction coefficient grid
    ExtinctionCoefficientGrid *extinction_grid = nullptr;
    if (do_extinction) {
      add_memory_allocation(
          ExtinctionCoefficientGrid::get_memory_size(number_of_angles),
          "ExtinctionCoefficientGrid", memory_log_file, memory_used);
      extinction_grid =
          new ExtinctionCoefficientGrid(number_of_angles, &theta[0]);
      quicksched.register_resource(*extinction_grid);
      quicksched.register_task(*extinction_grid);
      extinction_grid->link_resources(quicksched);
      tasks.push_back(extinction_grid);
    }
    //  - absorption coefficient grid
    AbsorptionCoefficientGrid *absorption_grid = nullptr;
    if (do_absorption) {
      add_memory_allocation(
          AbsorptionCoefficientGrid::get_memory_size(number_of_angles, ngauss),
          "AbsorptionCoefficientGrid", memory_log_file, memory_used);
      absorption_grid =
          new AbsorptionCoefficientGrid(number_of_angles, &theta[0], ngauss);
      quicksched.register_resource(*absorption_grid);
      quicksched.register_task(*absorption_grid);
      absorption_grid->link_resources(quicksched);
      tasks.push_back(absorption_grid);
    }
    //  - scattering matrix grid
    ScatteringMatrixGrid *scattering_grid = nullptr;
    if (do_scattering) {
      add_memory_allocation(ScatteringMatrixGrid::get_memory_size(
                                M_PI_2, number_of_angles, 2 * number_of_angles),
                            "ScatteringMatrixGrid", memory_log_file,
                            memory_used);
      scattering_grid = new ScatteringMatrixGrid(M_PI_2, number_of_angles,
                                                 2 * number_of_angles);
      quicksched.register_resource(*scattering_grid);
      quicksched.register_task(*scattering_grid);
      scattering_grid->link_resources(quicksched);
      tasks.push_back(scattering_grid);
    }

    //  - special Wigner D functions
    ExtinctionSpecialWignerDResources *extinction_special_wigner = nullptr;
    if (do_extinction) {
      add_memory_allocation(ExtinctionSpecialWignerDResources::get_memory_size(
                                _maximum_order, *extinction_grid),
                            "ExtinctionSpecialWignerDResources",
                            memory_log_file, memory_used);
      extinction_special_wigner = new ExtinctionSpecialWignerDResources(
          _maximum_order, *extinction_grid);
      quicksched.register_resource(*extinction_special_wigner);
      quicksched.register_task(*extinction_special_wigner);
      extinction_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*extinction_grid, *extinction_special_wigner);
      tasks.push_back(extinction_special_wigner);
    }

    AbsorptionSpecialWignerDResources *absorption_special_wigner = nullptr;
    if (do_absorption) {
      add_memory_allocation(AbsorptionSpecialWignerDResources::get_memory_size(
                                _maximum_order, *absorption_grid),
                            "AbsorptionSpecialWignerDResources",
                            memory_log_file, memory_used);
      absorption_special_wigner = new AbsorptionSpecialWignerDResources(
          _maximum_order, *absorption_grid);
      quicksched.register_resource(*absorption_special_wigner);
      quicksched.register_task(*absorption_special_wigner);
      absorption_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*absorption_grid, *absorption_special_wigner);
      tasks.push_back(absorption_special_wigner);
    }

    ScatteringMatrixSpecialWignerDResources *scattering_special_wigner =
        nullptr;
    if (do_scattering) {
      add_memory_allocation(
          ScatteringMatrixSpecialWignerDResources::get_memory_size(
              _maximum_order, *scattering_grid),
          "ScatteringMatrixSpecialWignerDResources", memory_log_file,
          memory_used);
      scattering_special_wigner = new ScatteringMatrixSpecialWignerDResources(
          _maximum_order, *scattering_grid);
      quicksched.register_resource(*scattering_special_wigner);
      quicksched.register_task(*scattering_special_wigner);
      scattering_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*scattering_grid, *scattering_special_wigner);
      tasks.push_back(scattering_special_wigner);
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
      add_memory_allocation(
          WignerDResources::get_memory_size(this_order, this_ngauss),
          "WignerDResources", memory_log_file, memory_used);
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
        add_memory_allocation(
            ParticleGeometryResource::get_memory_size(this_ngauss),
            "ParticleGeometryResource", memory_log_file, memory_used);
        ParticleGeometryResource *this_geometry =
            new ParticleGeometryResource(_shape_distribution.get_shape(is),
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
    if (do_extinction) {
      add_memory_allocation(
          total_number_of_Tmatrices *
              ExtinctionCoefficientResult::get_memory_size(number_of_angles),
          "ExtinctionCoefficientResult", memory_log_file, memory_used);
    }
    if (do_absorption) {
      add_memory_allocation(
          total_number_of_Tmatrices *
              AbsorptionCoefficientResult::get_memory_size(number_of_angles),
          "AbsorptionCoefficientResult", memory_log_file, memory_used);
    }
    if (do_scattering) {
      add_memory_allocation(total_number_of_Tmatrices *
                                ScatteringMatrixResult::get_memory_size(
                                    scattering_grid->get_number_of_angles()),
                            "ScatteringMatrixResult", memory_log_file,
                            memory_used);
    }

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
    add_memory_allocation(number_of_tmatrices * tmatrix_memory_requirement,
                          "TMatrixResources", memory_log_file, memory_used);
    if (verbose) {
      ctm_warning("Total memory usage: %s",
                  Utilities::human_readable_bytes(memory_used).c_str());
      ctm_warning("Space to store %" PRIuFAST32 " T-matrices.",
                  number_of_tmatrices);
    }

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
    std::vector<ExtinctionCoefficientResult *> unaveraged_extinction_results;
    if (do_extinction) {
      unaveraged_extinction_results.resize(total_number_of_Tmatrices, nullptr);
    }
    uint_fast32_t unaveraged_extinction_result_index = 0;
    std::vector<AbsorptionCoefficientResult *> unaveraged_absorption_results;
    if (do_absorption) {
      unaveraged_absorption_results.resize(total_number_of_Tmatrices, nullptr);
    }
    uint_fast32_t unaveraged_absorption_result_index = 0;
    std::vector<ScatteringMatrixResult *> unaveraged_scattering_results;
    if (do_scattering) {
      unaveraged_scattering_results.resize(total_number_of_Tmatrices, nullptr);
    }
    uint_fast32_t unaveraged_scattering_result_index = 0;
    const uint_fast32_t resource_offset = resources.size();
    resources.resize(resource_offset + total_number_of_interactions +
                         2 * total_number_of_Tmatrices,
                     nullptr);
    uint_fast32_t resource_index = resource_offset;
    results.resize(number_of_results * total_number_of_interactions, nullptr);
    uint_fast32_t result_index = 0;
    const uint_fast32_t task_offset = tasks.size();
    // for each T-matrix,
    // we need to add interaction and m=0 tasks for each quadrature point
    // we need to add m=/=0 tasks for each order
    // we need to add 1 alignment task and 3 reset tasks
    // we need to add an additional calculation task for each result
    const uint_fast32_t tasks_per_Tmatrix =
        2 * number_of_quadrature_tasks + _maximum_order + 4 + number_of_results;
    tasks.resize(task_offset +
                     number_of_results * total_number_of_interactions +
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
              _dust_properties.get_refractive_index(wavelength, particle_size,
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

          ExtinctionCoefficientResult *this_extinction_result = nullptr;
          ExtinctionShapeAveragingTask *extinction_averaging_task = nullptr;
          if (do_extinction) {
            this_extinction_result = new ExtinctionCoefficientResult(
                grain_type, particle_size, wavelength, number_of_angles);
            quicksched.register_resource(*this_extinction_result);
            results[result_index] = this_extinction_result;
            ++result_index;

            extinction_averaging_task = new ExtinctionShapeAveragingTask(
                _shape_distribution, *this_extinction_result);
            quicksched.register_task(*extinction_averaging_task);
            extinction_averaging_task->link_resources(quicksched);
            tasks[task_index] = extinction_averaging_task;
            ++task_index;
          }

          AbsorptionCoefficientResult *this_absorption_result = nullptr;
          AbsorptionShapeAveragingTask *absorption_averaging_task = nullptr;
          if (do_absorption) {
            this_absorption_result = new AbsorptionCoefficientResult(
                grain_type, particle_size, wavelength, number_of_angles);
            quicksched.register_resource(*this_absorption_result);
            results[result_index] = this_absorption_result;
            ++result_index;

            absorption_averaging_task = new AbsorptionShapeAveragingTask(
                _shape_distribution, *this_absorption_result);
            quicksched.register_task(*absorption_averaging_task);
            absorption_averaging_task->link_resources(quicksched);
            tasks[task_index] = absorption_averaging_task;
            ++task_index;
          }

          ScatteringMatrixResult *this_scattering_result = nullptr;
          ScatteringMatrixShapeAveragingTask *scattering_averaging_task =
              nullptr;
          if (do_scattering) {
            this_scattering_result = new ScatteringMatrixResult(
                grain_type, particle_size, wavelength,
                scattering_grid->get_number_of_angles());
            quicksched.register_resource(*this_scattering_result);
            results[result_index] = this_scattering_result;
            ++result_index;

            scattering_averaging_task = new ScatteringMatrixShapeAveragingTask(
                _shape_distribution, *this_scattering_result);
            quicksched.register_task(*scattering_averaging_task);
            scattering_averaging_task->link_resources(quicksched);
            tasks[task_index] = scattering_averaging_task;
            ++task_index;
          }

          // loop over all shapes
          for (uint_fast32_t ishape = 0; ishape < number_of_shapes; ++ishape) {

            ExtinctionCoefficientResult *this_unaveraged_extinction_result =
                nullptr;
            if (do_extinction) {
              this_unaveraged_extinction_result =
                  new ExtinctionCoefficientResult(grain_type, particle_size,
                                                  wavelength, number_of_angles);
              quicksched.register_resource(*this_unaveraged_extinction_result);
              unaveraged_extinction_results
                  [unaveraged_extinction_result_index] =
                      this_unaveraged_extinction_result;
              ++unaveraged_extinction_result_index;
              resources[resource_index] = this_unaveraged_extinction_result;
              ++resource_index;
            }

            AbsorptionCoefficientResult *this_unaveraged_absorption_result =
                nullptr;
            if (do_absorption) {
              this_unaveraged_absorption_result =
                  new AbsorptionCoefficientResult(grain_type, particle_size,
                                                  wavelength, number_of_angles);
              quicksched.register_resource(*this_unaveraged_absorption_result);
              unaveraged_absorption_results
                  [unaveraged_absorption_result_index] =
                      this_unaveraged_absorption_result;
              ++unaveraged_absorption_result_index;
              resources[resource_index] = this_unaveraged_absorption_result;
              ++resource_index;
            }

            ScatteringMatrixResult *this_unaveraged_scattering_result = nullptr;
            if (do_scattering) {
              this_unaveraged_scattering_result = new ScatteringMatrixResult(
                  grain_type, particle_size, wavelength,
                  scattering_grid->get_number_of_angles());
              quicksched.register_resource(*this_unaveraged_scattering_result);
              unaveraged_scattering_results
                  [unaveraged_scattering_result_index] =
                      this_unaveraged_scattering_result;
              ++unaveraged_scattering_result_index;
              resources[resource_index] = this_unaveraged_scattering_result;
              ++resource_index;
            }

            const OrientationDistribution &this_orientation =
                _alignment_distribution.get_distribution(
                    particle_size, _shape_distribution.get_shape(ishape));

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
                this_orientation, *this_single_Tmatrix, *this_ensemble_Tmatrix);
            quicksched.register_task(*alignment_task);
            alignment_task->link_resources(quicksched);
            tasks[task_index] = alignment_task;
            ++task_index;

            ExtinctionCoefficientTask *extinction_task = nullptr;
            if (do_extinction) {
              extinction_task = new ExtinctionCoefficientTask(
                  *extinction_grid, *this_interaction_variables,
                  *this_ensemble_Tmatrix, *nbased_resources,
                  *extinction_special_wigner,
                  *this_unaveraged_extinction_result);
              quicksched.register_task(*extinction_task);
              extinction_task->link_resources(quicksched);
              quicksched.link_tasks(*alignment_task, *extinction_task);
              quicksched.link_tasks(*extinction_special_wigner,
                                    *extinction_task);
              tasks[task_index] = extinction_task;
              ++task_index;

              extinction_averaging_task->add_input_coefficient(
                  quicksched, ishape, this_unaveraged_extinction_result);
              quicksched.link_tasks(*extinction_task,
                                    *extinction_averaging_task);
            }

            AbsorptionCoefficientTask *absorption_task = nullptr;
            if (do_absorption) {
              absorption_task = new AbsorptionCoefficientTask(
                  *absorption_grid, *this_interaction_variables,
                  *this_ensemble_Tmatrix, *nbased_resources,
                  *absorption_special_wigner,
                  *this_unaveraged_absorption_result, account_for_scattering);
              quicksched.register_task(*absorption_task);
              absorption_task->link_resources(quicksched);
              quicksched.link_tasks(*alignment_task, *absorption_task);
              quicksched.link_tasks(*absorption_special_wigner,
                                    *absorption_task);
              tasks[task_index] = absorption_task;
              ++task_index;

              absorption_averaging_task->add_input_coefficient(
                  quicksched, ishape, this_unaveraged_absorption_result);
              quicksched.link_tasks(*absorption_task,
                                    *absorption_averaging_task);
            }

            ScatteringMatrixTask *scattering_task = nullptr;
            if (do_scattering) {
              scattering_task = new ScatteringMatrixTask(
                  *scattering_grid, *this_interaction_variables,
                  *this_ensemble_Tmatrix, *nbased_resources,
                  *scattering_special_wigner,
                  *this_unaveraged_scattering_result);
              quicksched.register_task(*scattering_task);
              scattering_task->link_resources(quicksched);
              quicksched.link_tasks(*alignment_task, *scattering_task);
              quicksched.link_tasks(*scattering_special_wigner,
                                    *scattering_task);
              tasks[task_index] = scattering_task;
              ++task_index;

              scattering_averaging_task->add_input_matrices(
                  quicksched, ishape, this_unaveraged_scattering_result);
              quicksched.link_tasks(*scattering_task,
                                    *scattering_averaging_task);
            }

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
            if (extinction_task) {
              quicksched.link_tasks(*extinction_task, *reset_task2);
            }
            if (absorption_task) {
              quicksched.link_tasks(*absorption_task, *reset_task2);
            }
            if (scattering_task) {
              quicksched.link_tasks(*scattering_task, *reset_task2);
            }

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
