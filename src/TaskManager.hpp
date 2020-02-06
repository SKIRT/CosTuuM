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

  /**
   * @brief Class used to keep track of memory allocations.
   */
  class MemoryManager {
  private:
    /*! @brief Maximum allowed memory usage (in bytes). */
    const size_t _maximum_memory;

    /*! @brief Memory log file (if any). */
    std::ofstream *_memory_log_file;

    // memory management variables
    /*! @brief Memory already used (in bytes). */
    size_t _memory_used;

    // quicksched counter variables
    /*! @brief Number of QuickSched tasks. */
    size_t _number_of_tasks;

    /*! @brief Number of QuickSched resources. */
    size_t _number_of_resources;

    /*! @brief Number of QuickSched task dependencies. */
    size_t _number_of_task_dependencies;

    /*! @brief Number of QuickSched read/write resource dependencies. */
    size_t _number_of_resource_dependencies_readwrite;

    /*! @brief Number of QuickSched read only resource dependencies. */
    size_t _number_of_resource_dependencies_readonly;

    // pointer vector counter variables

    /*! @brief Number of Task pointers to store. */
    size_t _number_of_task_pointers;

    /*! @brief Number of Resource pointers to store. */
    size_t _number_of_resource_pointers;

    /*! @brief Number of Result pointers to store. */
    size_t _number_of_result_pointers;

    /*! @brief QuickSched memory usage. */
    size_t _quicksched_memory_usage;

  public:
    /**
     * @brief Constructor.
     *
     * @param maximum_memory Maximum allowed memory usage (in bytes).
     * @param write_memory_log Write a memory log file?
     * @param memory_log_file_name Name of the memory log file (if any).
     */
    inline MemoryManager(const size_t maximum_memory,
                         const bool write_memory_log = false,
                         const std::string memory_log_file_name = "")
        : _maximum_memory(maximum_memory), _memory_log_file(nullptr),
          _memory_used(0), _number_of_tasks(0), _number_of_resources(0),
          _number_of_task_dependencies(0),
          _number_of_resource_dependencies_readwrite(0),
          _number_of_resource_dependencies_readonly(0),
          _number_of_task_pointers(0), _number_of_resource_pointers(0),
          _number_of_result_pointers(0), _quicksched_memory_usage(0) {

      if (write_memory_log) {
        _memory_log_file = new std::ofstream(memory_log_file_name);
        *_memory_log_file << "# label\tallocation (bytes)\n";
        _memory_log_file->flush();
      }
    }

    /**
     * @brief Account for the given additional memory usage.
     *
     * @param size Additional memory usage size (in bytes).
     * @param label Label to mention in log file and error messages.
     */
    inline void add_memory_allocation(const size_t size,
                                      const std::string label) {

      if (_memory_log_file) {
        *_memory_log_file << label << "\t" << size << "\n";
        _memory_log_file->flush();
      }
      _memory_used += size;
      ctm_assert_message(_memory_used < _maximum_memory, "%s, %zu bytes",
                         label.c_str(), size);
    }

    /**
     * @brief Account for the given template Computable with the given
     * argument(s) and number of dependencies.
     *
     * @param label Label to show in the memory log and error messages.
     * @param number_of_dependent_tasks Number of QuickSched task dependencies
     * for this Computable.
     * @param args Additional arguments passed on to the
     * Computable::get_memory_size() function.
     * @tparam COMPUTABLE Computable class name.
     * @tparam ADDITIONAL_ARGUMENTS Variadic template to deal with arbitrary
     * additional parameters.
     */
    template <class COMPUTABLE, typename... ADDITIONAL_ARGUMENTS>
    inline void
    account_for_computable(const std::string label,
                           const uint_fast32_t number_of_dependent_tasks,
                           const ADDITIONAL_ARGUMENTS &... args) {

      add_memory_allocation(COMPUTABLE::get_memory_size(args...), label);
      ++_number_of_tasks;
      ++_number_of_resources;
      _number_of_resource_dependencies_readwrite +=
          COMPUTABLE::number_of_readwrite_resources();
      _number_of_resource_dependencies_readonly +=
          COMPUTABLE::number_of_readonly_resources();
      _number_of_task_dependencies += number_of_dependent_tasks;
      ++_number_of_task_pointers;
    }

    /**
     * @brief Account for the given template Task with the given argument(s) and
     * number of dependencies.
     *
     * @param label Label to show in the memory log and error messages.
     * @param number_of_dependent_tasks Number of QuickSched task dependencies
     * for this Task.
     * @param args Additional arguments passed on to the Task::get_memory_size()
     * function.
     * @tparam TASK Task class name.
     * @tparam ADDITIONAL_ARGUMENTS Variadic template to deal with arbitrary
     * additional parameters.
     */
    template <class TASK, typename... ADDITIONAL_ARGUMENTS>
    inline void account_for_task(const std::string label,
                                 const uint_fast32_t number_of_dependent_tasks,
                                 const ADDITIONAL_ARGUMENTS &... args) {

      add_memory_allocation(TASK::get_memory_size(args...), label);
      ++_number_of_tasks;
      _number_of_resource_dependencies_readwrite +=
          TASK::number_of_readwrite_resources();
      _number_of_resource_dependencies_readonly +=
          TASK::number_of_readonly_resources();
      _number_of_task_dependencies += number_of_dependent_tasks;
      ++_number_of_task_pointers;
    }

    /**
     * @brief Account for the given template Resource with the given
     * argument(s).
     *
     * @param label Label to show in the memory log and error messages.
     * @param args Additional arguments passed on to the
     * Resource::get_memory_size() function.
     * @tparam RESOURCE Resource class name.
     * @tparam ADDITIONAL_ARGUMENTS Variadic template to deal with arbitrary
     * additional parameters.
     */
    template <class RESOURCE, typename... ADDITIONAL_ARGUMENTS>
    inline void account_for_resource(const std::string label,
                                     const ADDITIONAL_ARGUMENTS &... args) {

      add_memory_allocation(RESOURCE::get_memory_size(args...), label);
      ++_number_of_resources;
      ++_number_of_resource_pointers;
    }

    /**
     * @brief Account for the given template Result with the given
     * argument(s).
     *
     * @param label Label to show in the memory log and error messages.
     * @param args Additional arguments passed on to the
     * Result::get_memory_size() function.
     * @tparam RESULT Result class name.
     * @tparam ADDITIONAL_ARGUMENTS Variadic template to deal with arbitrary
     * additional parameters.
     */
    template <class RESULT, typename... ADDITIONAL_ARGUMENTS>
    inline void account_for_result(const std::string label,
                                   const ADDITIONAL_ARGUMENTS &... args) {

      add_memory_allocation(RESULT::get_memory_size(args...), label);
      ++_number_of_resources;
      ++_number_of_result_pointers;
    }

    /**
     * @brief Get the number of task pointers.
     *
     * @return Number of task pointers.
     */
    inline size_t get_number_of_task_pointers() const {
      return _number_of_task_pointers;
    }

    /**
     * @brief Get the number of resource pointers.
     *
     * @return Number of resource pointers.
     */
    inline size_t get_number_of_resource_pointers() const {
      return _number_of_resource_pointers;
    }

    /**
     * @brief Get the number of result pointers.
     *
     * @return Number of result pointers.
     */
    inline size_t get_number_of_result_pointers() const {
      return _number_of_result_pointers;
    }

    /**
     * @brief Get the total memory usage.
     *
     * @return Total memory usage (in bytes).
     */
    inline size_t get_memory_used() const { return _memory_used; }

    /**
     * @brief Allocate sufficient memory to fit all QuickSched objects.
     *
     * @param quicksched QuickSched class wrapper.
     */
    inline void allocate_quicksched_memory(QuickSched &quicksched) {
      const size_t quicksched_memory = quicksched.reserve_memory(
          _number_of_tasks, _number_of_resources, _number_of_task_dependencies,
          _number_of_resource_dependencies_readwrite,
          _number_of_resource_dependencies_readonly);
      add_memory_allocation(quicksched_memory - _quicksched_memory_usage,
                            "QuickSched");
      _quicksched_memory_usage = quicksched_memory;
    }

    /**
     * @brief Display QuickSched object totals.
     */
    inline void output_totals() const {
      ctm_warning("Need to create %lu tasks, %lu resources, %lu task "
                  "dependencies, %lu read/write resource dependencies and %lu "
                  "read only resource dependencies",
                  _number_of_tasks, _number_of_resources,
                  _number_of_task_dependencies,
                  _number_of_resource_dependencies_readwrite,
                  _number_of_resource_dependencies_readonly);
    }
  };

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
    // number of quadrature and geometry tasks that we need
    const uint_fast32_t minimum_ngauss =
        _minimum_order * _gauss_legendre_factor;
    const uint_fast32_t maximum_ngauss =
        _maximum_order * _gauss_legendre_factor;
    const uint_fast32_t number_of_quadrature_tasks =
        _maximum_order - _minimum_order;
    const uint_fast32_t number_of_geometry_tasks =
        number_of_shapes * number_of_quadrature_tasks;

    // memory allocations: we first compute how much memory we will use
    // we need to keep track of memory to make sure we respect the user
    // defined limit
    // we will do a dummy creation of all tasks
    MemoryManager memory_manager(_maximum_memory_usage, write_memory_log,
                                 memory_log_file_name);
    // NBasedResources:
    //  task + resource
    memory_manager.account_for_computable<NBasedResources>("NBasedResources", 0,
                                                           _maximum_order);
    // GaussBasedResources:
    //  task + resource
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      memory_manager.account_for_computable<GaussBasedResources>(
          "GaussBasedResources", 0, this_ngauss);
    }
    // number of ExtinctionCoefficientGrids:
    //  task + resource
    if (do_extinction) {
      memory_manager.account_for_computable<ExtinctionCoefficientGrid>(
          "ExtinctionCoefficientGrid", 0, number_of_angles);
    }
    // number of AbsorptionCoefficientGrids:
    //  task + resource
    if (do_absorption) {
      memory_manager.account_for_computable<AbsorptionCoefficientGrid>(
          "AbsorptionCoefficientGrid", 0, number_of_angles, ngauss);
    }
    // number of ScatteringMatrixGrids:
    //  task + resource
    if (do_scattering) {
      memory_manager.account_for_computable<ScatteringMatrixGrid>(
          "ScatteringMatrixGrid", 0, M_PI_2, number_of_angles,
          2 * number_of_angles);
    }
    // number of ExtinctionSpecialWignerDResources:
    //  task + resource + 1 task dependency
    if (do_extinction) {
      memory_manager.account_for_computable<ExtinctionSpecialWignerDResources>(
          "ExtinctionSpecialWignerDResources", 1, _maximum_order,
          number_of_angles);
    }
    // number of AbsorptionSpecialWignerDResources:
    //  task + resource + 1 task dependency
    if (do_absorption) {
      memory_manager.account_for_computable<AbsorptionSpecialWignerDResources>(
          "AbsorptionSpecialWignerDResources", 1, _maximum_order,
          number_of_angles);
    }
    // number of ScatteringMatrixSpecialWignerDResources:
    //  task + resource + 1 task dependency
    if (do_scattering) {
      memory_manager
          .account_for_computable<ScatteringMatrixSpecialWignerDResources>(
              "ScatteringMatrixSpecialWignerDResources", 1, _maximum_order,
              number_of_angles);
    }
    // number of WignerDResources:
    //  task + resource + 1 task dependency
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_order = _minimum_order + i;
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      memory_manager.account_for_computable<WignerDResources>(
          "WignerDResources", 1, this_order, this_ngauss, this_order < 100);
    }
    // number of ParticleGeometryResources:
    //  task + resource + 1 task dependency
    for (uint_fast32_t is = 0; is < number_of_shapes; ++is) {
      for (uint_fast32_t ig = 0; ig < number_of_quadrature_tasks; ++ig) {
        const uint_fast32_t this_ngauss =
            minimum_ngauss + ig * _gauss_legendre_factor;
        memory_manager.account_for_computable<ParticleGeometryResource>(
            "ParticleGeometryResource", 1, this_ngauss);
      }
    }
    // at this point, we have accounted for all computables that need to be
    // allocated regardless of how much memory we have available
    // next, we do a dry run of the T-matrix task setup to account for the
    // resources that are required; this will give us the minimum required
    // memory usage to be able to execute any computation

    // loop over compositions
    for (uint_fast32_t icomp = 0; icomp < number_of_compositions; ++icomp) {
      // loop over all sizes
      for (uint_fast32_t isize = 0; isize < number_of_sizes; ++isize) {
        // loop over all wavelengths
        for (uint_fast32_t ilambda = 0; ilambda < number_of_wavelengths;
             ++ilambda) {
          memory_manager.account_for_resource<InteractionVariables>(
              "InteractionVariables");

          if (do_extinction) {
            memory_manager.account_for_result<ExtinctionCoefficientResult>(
                "ExtinctionCoefficientResult", number_of_angles);

            memory_manager.account_for_task<ExtinctionShapeAveragingTask>(
                "ExtinctionShapeAveragingTask", 0, _shape_distribution);
          }

          if (do_absorption) {
            memory_manager.account_for_result<AbsorptionCoefficientResult>(
                "AbsorptionCoefficientResult", number_of_angles);

            memory_manager.account_for_task<AbsorptionShapeAveragingTask>(
                "AbsorptionShapeAveragingTask", 0, _shape_distribution);
          }

          if (do_scattering) {
            memory_manager.account_for_result<ScatteringMatrixResult>(
                "ScatteringMatrixResult", number_of_angles);

            memory_manager.account_for_task<ScatteringMatrixShapeAveragingTask>(
                "ScatteringMatrixShapeAveragingTask", 0, _shape_distribution);
          }

          // loop over all shapes
          for (uint_fast32_t ishape = 0; ishape < number_of_shapes; ++ishape) {

            if (do_extinction) {
              memory_manager.account_for_resource<ExtinctionCoefficientResult>(
                  "ExtinctionCoefficientResult", number_of_angles);
            }

            if (do_absorption) {
              memory_manager.account_for_resource<AbsorptionCoefficientResult>(
                  "AbsorptionCoefficientResult", number_of_angles);
            }

            if (do_scattering) {
              memory_manager.account_for_resource<ScatteringMatrixResult>(
                  "ScatteringMatrixResult", number_of_angles);
            }

            // allocate the corresponding interaction variables and control
            // object
            memory_manager.account_for_resource<ConvergedSizeResources>(
                "ConvergedSizeResources");

            // loop over all orders
            for (uint_fast32_t ig = 0; ig < number_of_quadrature_tasks; ++ig) {
              // for safety, we overestimate the number of dependencies;
              // we don't know for sure whether this task will have an
              // additional dependency because it reuses a TMatrixResource
              memory_manager.account_for_task<InteractionTask>(
                  "InteractionTask", 2);

              // account for the m=0 task
              // again, we overestimate the number of dependencies
              memory_manager.account_for_task<TMatrixM0Task>("TMatrixM0Task",
                                                             4 + (ig > 0));
            }

            memory_manager.account_for_task<AlignmentAverageTask>(
                "AlignmentAverageTask", 0);

            if (do_extinction) {
              memory_manager.account_for_task<ExtinctionCoefficientTask>(
                  "ExtinctionCoefficientTask", 3);
            }

            if (do_absorption) {
              memory_manager.account_for_task<AbsorptionCoefficientTask>(
                  "AbsorptionCoefficientTask", 3);
            }

            if (do_scattering) {
              memory_manager.account_for_task<ScatteringMatrixTask>(
                  "ScatteringMatrixTask", 3);
            }

            memory_manager.account_for_task<ResetTMatrixResourceTask>(
                "ResetTMatrixResourceTask", 1);

            memory_manager.account_for_task<ResetTMatrixResourceTask>(
                "ResetTMatrixResourceTask",
                1 + do_extinction + do_absorption + do_scattering);

            memory_manager.account_for_task<ResetInteractionResourceTask>(
                "ResetInteractionResourceTask", 1);

            // loop over all m values to set up m=/=0 tasks
            for (uint_fast32_t i = 0; i < _maximum_order; ++i) {
              memory_manager.account_for_task<TMatrixMAllTask>(
                  "TMatrixMAllTask", 2);
            }
          }
        }
      }
    }

    // account for auxiliary space per thread
    memory_manager.add_memory_allocation(
        TMatrixAuxiliarySpaceManager::get_memory_size(
            quicksched.get_number_of_threads(), _maximum_order),
        "TMatrixAuxiliarySpaceManager");

    // account for the pointer storage
    memory_manager.add_memory_allocation(
        memory_manager.get_number_of_task_pointers() * sizeof(Task *), "Task*");
    memory_manager.add_memory_allocation(
        memory_manager.get_number_of_resource_pointers() * sizeof(Resource *),
        "Resource*");
    memory_manager.add_memory_allocation(
        memory_manager.get_number_of_result_pointers() * sizeof(Result *),
        "Result*");

    // we are done accounting for all necessary objects
    // first allocate QuickSched memory for all objects so far
    memory_manager.allocate_quicksched_memory(quicksched);

    // now allocate a single chunk of T-matrix and interaction resources to
    // get the per chunk memory requirements
    const size_t memory_size_required = memory_manager.get_memory_used();
    memory_manager.account_for_resource<TMatrixResource>("TMatrixResource",
                                                         _maximum_order);
    memory_manager.account_for_resource<TMatrixResource>("TMatrixResource",
                                                         _maximum_order);
    memory_manager.account_for_resource<InteractionResource>(
        "InteractionResource", _maximum_order, maximum_ngauss);
    // additional resource pointers
    memory_manager.add_memory_allocation(3 * sizeof(Resource *), "Resource*");
    memory_manager.allocate_quicksched_memory(quicksched);
    const size_t memory_per_Tmatrix =
        memory_manager.get_memory_used() - memory_size_required;
    // now figure out how many extra T-matrices we can store
    const size_t memory_left =
        _maximum_memory_usage - memory_manager.get_memory_used();
    const uint_fast32_t number_of_extra_tmatrices =
        std::min(memory_left / memory_per_Tmatrix, total_number_of_Tmatrices);
    for (uint_fast32_t itmatrix = 0; itmatrix < number_of_extra_tmatrices;
         ++itmatrix) {
      memory_manager.account_for_resource<TMatrixResource>("TMatrixResource",
                                                           _maximum_order);
      memory_manager.account_for_resource<TMatrixResource>("TMatrixResource",
                                                           _maximum_order);
      memory_manager.account_for_resource<InteractionResource>(
          "InteractionResource", _maximum_order, maximum_ngauss);
    }
    memory_manager.add_memory_allocation(3 * sizeof(Resource *), "Resource*");
    memory_manager.allocate_quicksched_memory(quicksched);

    const uint_fast32_t number_of_tmatrices = number_of_extra_tmatrices + 1;
    if (verbose) {
      ctm_warning(
          "Total memory usage: %s",
          Utilities::human_readable_bytes(memory_manager.get_memory_used())
              .c_str());
      ctm_warning("Space to store %" PRIuFAST32 " T-matrices.",
                  number_of_tmatrices);

      memory_manager.output_totals();
    }

    // now actually allocate and create everything

    // start with allocating the pointer vectors
    tasks.resize(memory_manager.get_number_of_task_pointers(), nullptr);
    resources.resize(memory_manager.get_number_of_resource_pointers(), nullptr);
    results.resize(memory_manager.get_number_of_result_pointers(), nullptr);
    // we will use running indices to fill them
    size_t running_task_index = 0;
    size_t running_resource_index = 0;
    size_t running_result_index = 0;

    // allocate auxiliary space per thread
    space_manager = new TMatrixAuxiliarySpaceManager(
        quicksched.get_number_of_threads(), _maximum_order);

    // step 1: set up the tasks that need to be done regardless of what the
    // parameter values are:

    //  - N based resources
    NBasedResources *nbased_resources = new NBasedResources(_maximum_order);
    tasks[running_task_index] = nbased_resources;
    ++running_task_index;
    quicksched.register_resource(*nbased_resources);
    quicksched.register_task(*nbased_resources);
    nbased_resources->link_resources(quicksched);

    //  - quadrature points
    std::vector<GaussBasedResources *> quadrature_points(
        number_of_quadrature_tasks, nullptr);
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      GaussBasedResources *this_quadrature_points =
          new GaussBasedResources(this_ngauss);
      quicksched.register_resource(*this_quadrature_points);
      quicksched.register_task(*this_quadrature_points);
      this_quadrature_points->link_resources(quicksched);
      quadrature_points[i] = this_quadrature_points;
      tasks[running_task_index] = this_quadrature_points;
      ++running_task_index;
    }

    //  - extinction coefficient grid
    ExtinctionCoefficientGrid *extinction_grid = nullptr;
    if (do_extinction) {
      extinction_grid =
          new ExtinctionCoefficientGrid(number_of_angles, &theta[0]);
      quicksched.register_resource(*extinction_grid);
      quicksched.register_task(*extinction_grid);
      extinction_grid->link_resources(quicksched);
      tasks[running_task_index] = extinction_grid;
      ++running_task_index;
    }
    //  - absorption coefficient grid
    AbsorptionCoefficientGrid *absorption_grid = nullptr;
    if (do_absorption) {
      absorption_grid =
          new AbsorptionCoefficientGrid(number_of_angles, &theta[0], ngauss);
      quicksched.register_resource(*absorption_grid);
      quicksched.register_task(*absorption_grid);
      absorption_grid->link_resources(quicksched);
      tasks[running_task_index] = absorption_grid;
      ++running_task_index;
    }
    //  - scattering matrix grid
    ScatteringMatrixGrid *scattering_grid = nullptr;
    if (do_scattering) {
      scattering_grid = new ScatteringMatrixGrid(M_PI_2, number_of_angles,
                                                 2 * number_of_angles);
      quicksched.register_resource(*scattering_grid);
      quicksched.register_task(*scattering_grid);
      scattering_grid->link_resources(quicksched);
      tasks[running_task_index] = scattering_grid;
      ++running_task_index;
    }

    //  - special Wigner D functions
    ExtinctionSpecialWignerDResources *extinction_special_wigner = nullptr;
    if (do_extinction) {
      extinction_special_wigner = new ExtinctionSpecialWignerDResources(
          _maximum_order, *extinction_grid);
      quicksched.register_resource(*extinction_special_wigner);
      quicksched.register_task(*extinction_special_wigner);
      extinction_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*extinction_grid, *extinction_special_wigner);
      tasks[running_task_index] = extinction_special_wigner;
      ++running_task_index;
    }

    AbsorptionSpecialWignerDResources *absorption_special_wigner = nullptr;
    if (do_absorption) {
      absorption_special_wigner = new AbsorptionSpecialWignerDResources(
          _maximum_order, *absorption_grid);
      quicksched.register_resource(*absorption_special_wigner);
      quicksched.register_task(*absorption_special_wigner);
      absorption_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*absorption_grid, *absorption_special_wigner);
      tasks[running_task_index] = absorption_special_wigner;
      ++running_task_index;
    }

    ScatteringMatrixSpecialWignerDResources *scattering_special_wigner =
        nullptr;
    if (do_scattering) {
      scattering_special_wigner = new ScatteringMatrixSpecialWignerDResources(
          _maximum_order, *scattering_grid);
      quicksched.register_resource(*scattering_special_wigner);
      quicksched.register_task(*scattering_special_wigner);
      scattering_special_wigner->link_resources(quicksched);
      quicksched.link_tasks(*scattering_grid, *scattering_special_wigner);
      tasks[running_task_index] = scattering_special_wigner;
      ++running_task_index;
    }

    //  - Wigner D functions
    std::vector<WignerDResources *> wignerdm0(number_of_quadrature_tasks,
                                              nullptr);
    for (uint_fast32_t i = 0; i < number_of_quadrature_tasks; ++i) {
      const uint_fast32_t this_order = _minimum_order + i;
      const uint_fast32_t this_ngauss =
          minimum_ngauss + i * _gauss_legendre_factor;
      WignerDResources *this_wignerdm0 = new WignerDResources(
          this_order, this_ngauss, *quadrature_points[i], this_order < 100);
      quicksched.register_resource(*this_wignerdm0);
      quicksched.register_task(*this_wignerdm0);
      this_wignerdm0->link_resources(quicksched);
      quicksched.link_tasks(*quadrature_points[i], *this_wignerdm0);
      wignerdm0[i] = this_wignerdm0;
      tasks[running_task_index] = this_wignerdm0;
      ++running_task_index;
    }

    // step 2: compute shape quadrature points and generate shape based
    // resources that are shared between all parameter values
    std::vector<ParticleGeometryResource *> geometries(number_of_geometry_tasks,
                                                       nullptr);
    for (uint_fast32_t is = 0; is < number_of_shapes; ++is) {
      for (uint_fast32_t ig = 0; ig < number_of_quadrature_tasks; ++ig) {
        const uint_fast32_t this_ngauss =
            minimum_ngauss + ig * _gauss_legendre_factor;
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
        tasks[running_task_index] = this_geometry;
        ++running_task_index;
      }
    }

    // step 4: loop over all parameter values and set up parameter specific
    // tasks

    // now actually allocate the resources
    std::vector<InteractionResource *> interaction_resources(
        number_of_tmatrices, nullptr);
    std::vector<TMatrixResource *> tmatrices(2 * number_of_tmatrices, nullptr);
    for (uint_fast32_t i = 0; i < number_of_tmatrices; ++i) {
      interaction_resources[i] =
          new InteractionResource(_maximum_order, maximum_ngauss);
      quicksched.register_resource(*interaction_resources[i]);
      resources[running_resource_index] = interaction_resources[i];
      ++running_resource_index;

      tmatrices[2 * i] = new TMatrixResource(_maximum_order);
      resources[running_resource_index] = tmatrices[2 * i];
      ++running_resource_index;
      quicksched.register_resource(*tmatrices[2 * i]);
      for (uint_fast32_t m = 0; m < _maximum_order + 1; ++m) {
        quicksched.register_resource(tmatrices[2 * i]->get_m_resource(m),
                                     tmatrices[2 * i]);
      }
      tmatrices[2 * i + 1] = new TMatrixResource(_maximum_order);
      resources[running_resource_index] = tmatrices[2 * i + 1];
      ++running_resource_index;
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
    uint_fast32_t resource_reuse_index = 0;

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
          resources[running_resource_index] = this_interaction_variables;
          ++running_resource_index;

          ExtinctionCoefficientResult *this_extinction_result = nullptr;
          ExtinctionShapeAveragingTask *extinction_averaging_task = nullptr;
          if (do_extinction) {
            this_extinction_result = new ExtinctionCoefficientResult(
                grain_type, particle_size, wavelength, number_of_angles);
            quicksched.register_resource(*this_extinction_result);
            results[running_result_index] = this_extinction_result;
            ++running_result_index;

            extinction_averaging_task = new ExtinctionShapeAveragingTask(
                _shape_distribution, *this_extinction_result);
            quicksched.register_task(*extinction_averaging_task);
            extinction_averaging_task->link_resources(quicksched);
            tasks[running_task_index] = extinction_averaging_task;
            ++running_task_index;
          }

          AbsorptionCoefficientResult *this_absorption_result = nullptr;
          AbsorptionShapeAveragingTask *absorption_averaging_task = nullptr;
          if (do_absorption) {
            this_absorption_result = new AbsorptionCoefficientResult(
                grain_type, particle_size, wavelength, number_of_angles);
            quicksched.register_resource(*this_absorption_result);
            results[running_result_index] = this_absorption_result;
            ++running_result_index;

            absorption_averaging_task = new AbsorptionShapeAveragingTask(
                _shape_distribution, *this_absorption_result);
            quicksched.register_task(*absorption_averaging_task);
            absorption_averaging_task->link_resources(quicksched);
            tasks[running_task_index] = absorption_averaging_task;
            ++running_task_index;
          }

          ScatteringMatrixResult *this_scattering_result = nullptr;
          ScatteringMatrixShapeAveragingTask *scattering_averaging_task =
              nullptr;
          if (do_scattering) {
            this_scattering_result = new ScatteringMatrixResult(
                grain_type, particle_size, wavelength,
                scattering_grid->get_number_of_angles());
            quicksched.register_resource(*this_scattering_result);
            results[running_result_index] = this_scattering_result;
            ++running_result_index;

            scattering_averaging_task = new ScatteringMatrixShapeAveragingTask(
                _shape_distribution, *this_scattering_result);
            quicksched.register_task(*scattering_averaging_task);
            scattering_averaging_task->link_resources(quicksched);
            tasks[running_task_index] = scattering_averaging_task;
            ++running_task_index;
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
              resources[running_resource_index] =
                  this_unaveraged_extinction_result;
              ++running_resource_index;
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
              resources[running_resource_index] =
                  this_unaveraged_absorption_result;
              ++running_resource_index;
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
              resources[running_resource_index] =
                  this_unaveraged_scattering_result;
              ++running_resource_index;
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
            resources[running_resource_index] = this_converged_size;
            ++running_resource_index;

            InteractionResource *this_interaction_resource =
                interaction_resources[resource_reuse_index];
            TMatrixResource *this_single_Tmatrix =
                tmatrices[2 * resource_reuse_index];
            TMatrixResource *this_ensemble_Tmatrix =
                tmatrices[2 * resource_reuse_index + 1];

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
              if (unlock_tasks[resource_reuse_index] != nullptr) {
                quicksched.link_tasks(*unlock_tasks[resource_reuse_index],
                                      *this_interaction);
              }
              tasks[running_task_index] = this_interaction;
              ++running_task_index;

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
                if (unlock_tasks[resource_reuse_index] != nullptr) {
                  quicksched.link_tasks(*unlock_tasks[resource_reuse_index],
                                        *this_m0task);
                }
              }
              tasks[running_task_index] = this_m0task;
              ++running_task_index;
              previous_m0task = this_m0task;
            }

            AlignmentAverageTask *alignment_task = new AlignmentAverageTask(
                this_orientation, *this_single_Tmatrix, *this_ensemble_Tmatrix);
            quicksched.register_task(*alignment_task);
            alignment_task->link_resources(quicksched);
            tasks[running_task_index] = alignment_task;
            ++running_task_index;

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
              tasks[running_task_index] = extinction_task;
              ++running_task_index;

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
              tasks[running_task_index] = absorption_task;
              ++running_task_index;

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
              tasks[running_task_index] = scattering_task;
              ++running_task_index;

              scattering_averaging_task->add_input_matrices(
                  quicksched, ishape, this_unaveraged_scattering_result);
              quicksched.link_tasks(*scattering_task,
                                    *scattering_averaging_task);
            }

            ResetTMatrixResourceTask *reset_task1 =
                new ResetTMatrixResourceTask(*this_single_Tmatrix);
            quicksched.register_task(*reset_task1);
            reset_task1->link_resources(quicksched);
            tasks[running_task_index] = reset_task1;
            ++running_task_index;
            quicksched.link_tasks(*alignment_task, *reset_task1);

            ResetTMatrixResourceTask *reset_task2 =
                new ResetTMatrixResourceTask(*this_ensemble_Tmatrix);
            quicksched.register_task(*reset_task2);
            reset_task2->link_resources(quicksched);
            tasks[running_task_index] = reset_task2;
            ++running_task_index;
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
            tasks[running_task_index] = reset_task3;
            ++running_task_index;
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
              tasks[running_task_index] = this_malltask;
              ++running_task_index;
            }

            unlock_tasks[resource_reuse_index] = reset_task3;
            ++resource_reuse_index;
            resource_reuse_index %= number_of_tmatrices;
          }
        }
      }
    }
  }
};

#endif // TASKMANAGER_HPP
