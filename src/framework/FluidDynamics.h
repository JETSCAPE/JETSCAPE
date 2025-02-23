/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// This is a general basic class for hydrodynamics

#ifndef FLUIDDYNAMICS_H
#define FLUIDDYNAMICS_H

#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include "FluidCellInfo.h"
#include "FluidEvolutionHistory.h"
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "LiquefierBase.h"
#include "PreequilibriumDynamics.h"
#include "RealType.h"
#include "SurfaceCellInfo.h"

namespace Jetscape {

/**
 * @brief Flags for hydrodynamics status.
 *
 * This enum defines the various states that the hydrodynamics process can be
 * in. It helps track the progress and status of the simulation or computation.
 *
 * @enum HydroStatus
 */
enum HydroStatus { NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR };

/**
 * @class Parameter
 * @brief A helper class for managing hydro parameters file name.
 *
 * This class is designed to hold and manage the file name of the hydro input
 * parameters. It contains a single member variable that stores the file name as
 * a string.
 */
class Parameter {
 public:
  /**
   * @brief The file name for the hydro input parameters.
   */
  char *hydro_input_filename;
};

/**
 * @class FluidDynamics
 * @brief Represents a module for hydrodynamic evolution in high-energy nuclear
 * collisions.
 *
 * This class is derived from `JetScapeModuleBase` and provides functionality
 * for initializing and evolving hydrodynamic fields, managing evolution
 * history, and querying fluid properties such as energy density, entropy
 * density, temperature, and quark-gluon plasma fraction. It also supports
 * interfacing with other modules through a polymorphic structure and provides
 * customization hooks for different hydrodynamic models.
 *
 * The class includes:
 * - Parameters for initialization and freeze-out conditions
 * - Evolution history tracking and surface cell management
 * - Accessor methods for fluid properties and state
 * - Integration with external modules for liquefaction and jet energy
 * deposition
 *
 * @see JetScapeModuleBase, InitialState, PreequilibriumDynamics, LiquefierBase,
 * EvolutionHistory
 */
class FluidDynamics : public JetScapeModuleBase {
 protected:
  /**
   * @brief Record the start and end proper time of the hydro simulation.
   *
   * These variables represent the start and end proper time in (fm/c).
   *
   * @var Jetscape::real hydro_tau_0
   * @var Jetscape::real hydro_tau_max
   */
  Jetscape::real hydro_tau_0, hydro_tau_max;

  /**
   * @brief Record the hydro freeze out temperature.
   *
   * This variable stores the hydro freeze out temperature in GeV.
   *
   * @var Jetscape::real hydro_freeze_out_temperature
   */
  Jetscape::real hydro_freeze_out_temperature;

  /**
   * @brief Record the status of the hydrodynamics.
   *
   * @var HydroStatus hydro_status
   */
  HydroStatus hydro_status;

  // add initial state shared pointer
  /** A shared pointer to the initial state class. */
  std::shared_ptr<InitialState> ini;
  /** A shared pointer to the pre-equilibrium dynamics class. */
  std::shared_ptr<PreequilibriumDynamics> pre_eq_ptr;

  double eta;
  bool boost_invariant_;
  Parameter parameter_list;

  // How to store this data? In memory or hard disk?
  // 3D hydro may eat out the memory,
  // for large dataset, std::deque is better than std::vector.

  /** Stores the evolution history. */
  EvolutionHistory bulk_info;

  /** Stores the surface cell information. */
  std::vector<SurfaceCellInfo> surfaceCellVector_;

  /** A weak pointer to the liquefier object. */
  std::weak_ptr<LiquefierBase> liquefier_ptr;

 public:
  /**
   * @brief Default constructor for the FluidDynamics class.
   *
   * Initializes the object with a default task ID of "FluidDynamics" and
   * sets the eta parameter to -99.99. Additionally, the boost_invariant_
   * is set to true. This constructor is typically used to initialize
   * an instance of FluidDynamics with default values.
   *
   * @note The task ID is set via the SetId method to "FluidDynamics".
   * @note The eta parameter is initialized to -99.99 to indicate an
   * uninitialized state.
   */
  FluidDynamics();

  /**
   * @brief Default destructor for the FluidDynamics class.
   *
   * This destructor is responsible for cleaning up any resources held by
   * the FluidDynamics object before it is destroyed. It performs necessary
   * disconnections and logs verbose output to indicate the destruction process.
   *
   * @note This destructor calls `disconnect_all()` to ensure that all
   *       connections are properly severed when the object is destroyed.
   *
   * @see disconnect_all()
   */
  virtual ~FluidDynamics();

  /**
   * @brief Initializes the FluidDynamics module by reading input parameters
   *        from the XML file under the <Hydro> tag, retrieving the initial
   * state physics information, and initializing related modules and tasks.
   *
   * This function uses the `JetScapeSignalManager` instance to retrieve the
   * initial state physics data and pre-equilibrium module. It then calls the
   * `InitializeHydro(parameter_list)` function and the `InitTask()` function
   * to initialize various components within the current module. This explicit
   * call is useful for initializing modules like @a Brick, @a MpiMusic, or
   * @a OSU-HYDRO if attached as a @a polymorphic class. Additionally, the
   * function initializes tasks within the current module using
   * `JetScapeTask::InitTasks()`.
   *
   * @note This function assumes the availability of the `parameter_list` and
   *       initializes the module by using its specific tasks.
   *
   * @sa For more information on polymorphism in C++, refer to C++ polymorphism
   * documentation.
   *
   * @see JetScapeSignalManager, JetScapeModuleBase, InitializeHydro, InitTask
   */
  virtual void Init();

  /**
   * @brief Executes the hydrodynamic evolution and additional tasks within the
   * current module.
   *
   * This function explicitly calls the `EvolveHydro()` method to perform the
   * hydrodynamic evolution, which is defined in various modules like @a Brick,
   * @a MpiMusic, or @a OSU-HYDRO if they are attached as polymorphic classes.
   * In addition, the function also executes the tasks associated with the
   * current module.
   *
   * It provides verbose logging for debugging purposes, including information
   * about the system state such as entropy density distribution size and event
   * number.
   *
   * @sa Read about @a polymorphism in C++ for more information on how
   * polymorphic classes are used in the context of hydrodynamic evolution.
   *
   * @note This function is designed to be virtual, allowing for potential
   * overrides in derived classes to implement specific behavior.
   *
   * @see EvolveHydro()
   * @see JetScapeTask::ExecuteTasks()
   */
  virtual void Exec();

  /**
   * @brief Clears the fluid dynamics data and resets related components.
   *
   * This function performs the following operations:
   * - Clears the evolution data by calling `clear_up_evolution_data()`.
   * - If the `liquefier_ptr` is not uninitialized, it locks the `weak_ptr` and
   * calls `Clear()` on the `liquefier` object to clear its data as well.
   *
   * @note This function is virtual and may be overridden in derived classes.
   */
  virtual void Clear();

  /**
   * @brief Retrieves the parameter list for the fluid dynamics task.
   *
   * This function returns a reference to the `Parameter` object, which contains
   * a file name for the fluid dynamics task.
   *
   * @return Parameter& A reference to the `Parameter` object.
   *
   * @sa Implementation of the class Parameter.
   */
  Parameter &GetParameterList() { return (parameter_list); }

  /**
   * @brief Collects header information for writer modules.
   *
   * This function extracts and sets the event plane angle in the header of the
   * provided JetScapeWriter instance, if the writer is valid.
   *
   * @param w A weak pointer to a JetScapeWriter instance. This pointer is used
   * to retrieve the header for setting the event plane angle.
   *
   * @note If the weak pointer is expired or invalid, no operation is performed.
   */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Retrieves the generated event plane angle.
   *
   * This function returns a placeholder value of -999 by default and is
   * expected to be overridden by derived classes that provide a meaningful
   * implementation for retrieving the event plane angle.
   *
   * @return The event plane angle as a double. Default is -999 if not
   * overridden.
   */
  virtual double GetEventPlaneAngle() { return (-999); }

  /**
   * @brief Stores the temperature of the fluid cell at a given location.
   *
   * This function retrieves the temperature of the fluid cell at the specified
   * coordinates (t or tau, x, y, z or eta) and stores it in the input variable
   * @p mT. It can be overridden by modules attached to the FluidDynamics class.
   *
   * @param t  Time or proper time (tau) coordinate.
   * @param x  Spatial x-coordinate.
   * @param y  Spatial y-coordinate.
   * @param z  Spatial z-coordinate or rapidity (eta).
   * @param[out] mT Variable to store the temperature.
   */
  virtual void GetTemperature(double t, double x, double y, double z,
                              double &mT) {
    mT = GetTemperature(t, x, y, z);
  }

  /**
   * @brief Retrieves the properties of the fluid cell at a given location.
   *
   * This function calls GetHydroInfo(t, x, y, z, fCell) to obtain the
   * properties of the fluid cell at the specified coordinates. It can be
   * overridden by modules attached to the FluidDynamics class.
   *
   * @param t     The time or proper time (tau) coordinate.
   * @param x     The spatial x coordinate.
   * @param y     The spatial y coordinate.
   * @param z     The spatial z coordinate or rapidity (eta).
   * @param fCell A unique pointer to a FluidCellInfo object where the retrieved
   *              data will be stored.
   */
  virtual void GetHydroCell(double t, double x, double y, double z,
                            std::unique_ptr<FluidCellInfo> &fCell) {
    GetHydroInfo(t, x, y, z, fCell);
  }

  // currently we have no standard for passing configurations
  // pure virtual function; to be implemented by users
  // should make it easy to save evolution history to bulk_info

  /**
   * @brief Default function to initialize the hydrodynamics.
   *
   * This pure virtual function serves as a default initialization method for
   * hydrodynamics and can be overridden by different modules as needed.
   * Implementing this function should facilitate saving evolution history
   * to bulk_info.
   *
   * @note Currently, there is no standard for passing configurations.
   *
   * @param parameter_list An object of the class Parameter containing
   *        initialization parameters.
   */
  virtual void InitializeHydro(Parameter parameter_list){};

  /**
   * @brief Default function to evolve the hydrodynamics.
   *
   * This function represents the default implementation for evolving
   * the hydrodynamics. It can be overridden by different modules.
   *
   * @note This is a virtual function, allowing for customization in
   *       derived classes.
   */
  virtual void EvolveHydro(){};

  /**
   * @brief Retrieves the current status of the hydrodynamics.
   *
   * This function returns the status of the hydrodynamics system, indicating
   * whether it has not started, is initialized, is evolving, has finished, or
   * encountered an error.
   *
   * @return int The status of the hydrodynamics:
   *   - NOT_START: The hydrodynamics process has not started.
   *   - INITIALIZED: The hydrodynamics process is initialized.
   *   - EVOLVING: The hydrodynamics process is in progress.
   *   - FINISHED: The hydrodynamics process is completed.
   *   - ERROR: An error occurred during the hydrodynamics process.
   */
  int GetHydroStatus() const { return (hydro_status); }

  /**
   * @brief Stores the fluid cell evolution history in the bulk information.
   *
   * This function takes a unique pointer to a `FluidCellInfo` object and stores
   * its content in the `data` member of the `bulk_info` object. The function
   * dereferences the unique pointer and pushes the content to the `data`
   * vector.
   *
   * @param[in,out] fluid_cell_info_ptr A unique pointer to the `FluidCellInfo`
   * object that contains the data to be stored. The object will not be
   * modified, but the pointer will not be used after the function call.
   */
  void StoreHydroEvolutionHistory(
      std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
    bulk_info.data.push_back(*fluid_cell_info_ptr);
  }

  /**
   * @brief Stores a surface cell information in the surface cell vector.
   *
   * This function adds the provided `SurfaceCellInfo` object to the
   * `surfaceCellVector_`. It allows the surface cell information to be stored
   * for later use or processing.
   *
   * @param[in] surface_cell_info The surface cell information to be stored.
   * The object should contain all necessary details of the surface cell that
   * need to be recorded.
   *
   * @note This function modifies the internal `surfaceCellVector_` by pushing
   * the provided `surface_cell_info` to the back of the vector.
   */
  void StoreSurfaceCell(SurfaceCellInfo &surface_cell_info) {
    surfaceCellVector_.push_back(surface_cell_info);
  }

  /**
   * @brief Retrieves the surface cell vector and logs its size.
   *
   * This function copies the contents of the internal `surfaceCellVector_` to
   * the provided `surfCellVec` vector. After copying, it logs the size of the
   * vector with a message.
   *
   * @param[in, out] surfCellVec The vector to which the surface cell
   * information is copied. It will contain the surface cells after the function
   * call.
   *
   * @note This function does not modify the `surfaceCellVector_` but rather
   * copies its contents to the provided `surfCellVec`.
   */
  void getSurfaceCellVector(std::vector<SurfaceCellInfo> &surfCellVec) {
    surfCellVec = surfaceCellVector_;
    JSINFO << "Fluid out: surface vector size = " << surfCellVec.size();
  }

  /**
   * @brief Clears the surface cell vector.
   *
   * This function clears the contents of the `surfaceCellVector_`, effectively
   * removing all elements from the vector.
   *
   * @note The operation does not return any value and operates directly on the
   *       `surfaceCellVector_` member variable.
   *
   * @see surfaceCellVector_
   */
  void clearSurfaceCellVector() { surfaceCellVector_.clear(); }

  /**
   * @brief Clears up evolution data from the bulk information.
   *
   * This function invokes the `clear_up_evolution_data` method of the
   * `bulk_info` object to clear or reset any evolution-related data it holds.
   * The exact implementation details of how the data is cleared are dependent
   * on the `bulk_info` class and its `clear_up_evolution_data` method.
   *
   * @note This function assumes that `bulk_info` is an object of a class that
   * contains a valid `clear_up_evolution_data` method.
   */
  void clear_up_evolution_data() { bulk_info.clear_up_evolution_data(); }

  /**
   * @brief Retrieves the start time (or tau) for hydrodynamic evolution.
   *
   * This function sets the provided reference parameter `tau0` to the start
   * time for the hydrodynamic evolution, which is represented by the member
   * variable `hydro_tau_0`.
   *
   * @param[out] tau0 A reference to a double where the start time (tau) will be
   * stored.
   *
   * @return void
   */
  void GetHydroStartTime(double &tau0) { tau0 = hydro_tau_0; }

  /**
   * @brief Returns the end time (or tau) for the hydrodynamic evolution.
   *
   * This function provides the maximum value of the hydrodynamic evolution
   * time, represented as the end time (or tau), for the simulation or
   * calculation.
   *
   * @return The end time (or tau) for the hydrodynamic evolution.
   */
  Jetscape::real GetHydroEndTime() const { return (hydro_tau_max); }

  /**
   * @brief Returns the freeze-out temperature.
   *
   * This function retrieves the current freeze-out temperature value,
   * which is typically used in hydrodynamic simulations.
   *
   * @return Freeze-out temperature (of type Jetscape::real).
   */
  Jetscape::real GetHydroFreezeOutTemperature() const {
    return (hydro_freeze_out_temperature);
  }

  /**
   * @brief Retrieves the hydro information at a given space-time point.
   *
   * This function fetches the hydro information at the specified space-time
   * point (t, x, y, z or eta). It checks if the time and spatial coordinates
   * are within the range of the evolution history. If they are out of range, it
   * throws an `InvalidSpaceTimeRange` exception.
   *
   * @param t The time or tau coordinate.
   * @param x The spatial coordinate (x).
   * @param y The spatial coordinate (y).
   * @param z The spatial or eta coordinate (z or eta).
   * @param fluid_cell_info_ptr A unique pointer to the `FluidCellInfo` class
   * that will store the retrieved hydro information.
   *
   * @throws std::runtime_error If hydro evolution is not finished or if the
   * evolution history is empty.
   * @throws InvalidSpaceTimeRange If the provided space-time coordinates are
   * out of the evolution history range.
   */
  virtual void GetHydroInfo(
      Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
      std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
    if (hydro_status != FINISHED || bulk_info.data.size() == 0) {
      throw std::runtime_error(
          "Hydro evolution is not finished "
          "or EvolutionHistory is empty");
    }
    // judge whether to use 2D interpolation or 3D interpolation
    if (!bulk_info.tau_eta_is_tz) {
      Jetscape::real tau = std::sqrt(t * t - z * z);
      Jetscape::real eta = 0.5 * (std::log(t + z) - std::log(t - z));
      bulk_info.CheckInRange(tau, x, y, eta);
      // return bulk_info.get(tau, x, y, eta);
    } else {
      bulk_info.CheckInRange(t, x, y, z);
      // return bulk_info.get(t, x, y, z);
    }
  }

  /**
   * @brief Prints out the information of the fluid cell to the screen.
   *
   * This function outputs detailed information about the fluid cell's state,
   * including various physical properties such as energy density, entropy
   * density, temperature, pressure, chemical potentials, velocity components,
   * and shear/bulk viscosities.
   *
   * @param fluid_cell_info_ptr A pointer to the `FluidCellInfo` object
   * containing the fluid cell's data.
   *
   * @note This function outputs the information using the `JSINFO` logging
   * system. The values printed are in appropriate units such as GeV/fm^3 for
   * energy and pressure, and GeV for chemical potentials.
   *
   * @see FluidCellInfo
   */
  void PrintFluidCellInformation(FluidCellInfo *fluid_cell_info_ptr);

  /**
   * @brief Calculates the hypersurface for the Cooper-Frye or recombination
   * model.
   *
   * This function is designed to return the hypersurface at a given constant
   * temperature and can be overridden by different modules for more specific
   * implementations. It is part of the fluid dynamics computations, where the
   * hypersurface is needed to track the evolution of the system at the
   * specified temperature.
   *
   * @param T_sw The switching temperature at which the hypersurface is
   * computed. This is a scalar value representing the temperature threshold.
   * @param surface_cells A reference to a vector of `SurfaceCellInfo` objects
   * that will be populated with the resulting surface cells. These cells
   * represent the hypersurface computed at the given temperature.
   *
   * @return void This function does not return any value. It modifies the
   * `surface_cells` vector with the calculated hypersurface information.
   *
   * @note The detailed implementation of how the hypersurface is calculated is
   * left to the hydro developer. The default behavior calls the `SurfaceFinder`
   * to determine the surface and populate the `surface_cells`.
   *
   * @see SurfaceFinder
   * @see SurfaceCellInfo
   */
  void FindAConstantTemperatureSurface(
      Jetscape::real T_sw, std::vector<SurfaceCellInfo> &surface_cells);

  // all the following functions will call function GetHydroInfo()
  // to get thermaldynamic and dynamical information at a space-time point
  // (time, x, y, z)

  /**
   * @brief Returns the energy density at a given space-time point.
   *
   * This function retrieves the thermodynamic and dynamical information at a
   * specified space-time point (time, x, y, z or eta) by calling the
   * `GetHydroInfo()` function, which provides the relevant fluid cell data. It
   * then extracts and returns the energy density [GeV] from that information.
   *
   * @param time The time or tau coordinate (in physical units of time or
   * rapidity).
   * @param x The space coordinate (in physical units of length).
   * @param y The space coordinate (in physical units of length).
   * @param z The space or eta coordinate (in physical units of length or
   * rapidity).
   *
   * @return The energy density at the given space-time point, in GeV.
   *
   * @see GetHydroInfo()
   */
  virtual Jetscape::real GetEnergyDensity(Jetscape::real time, Jetscape::real x,
                                          Jetscape::real y, Jetscape::real z);

  /**
   * @brief Returns the entropy density at a given space-time point.
   *
   * This function computes and returns the entropy density at a specified
   * space-time point (time, x, y, z). It retrieves the fluid cell information
   * from the given space-time coordinates and accesses the entropy density
   * value stored within the fluid cell.
   *
   * @param time The time or tau coordinate of the space-time point.
   * @param x The x-space coordinate of the space-time point.
   * @param y The y-space coordinate of the space-time point.
   * @param z The z-space coordinate or eta coordinate of the space-time point.
   *
   * @return The entropy density at the given space-time point, in units of
   * [GeV].
   *
   * @see GetHydroInfo()
   */
  virtual Jetscape::real GetEntropyDensity(Jetscape::real time,
                                           Jetscape::real x, Jetscape::real y,
                                           Jetscape::real z);

  /**
   * @brief Returns the temperature at a specific spacetime point.
   *
   * This function calculates and returns the temperature at the given spacetime
   * coordinates (time, x, y, z) in GeV. It retrieves the necessary fluid
   * dynamics information and extracts the temperature at the specified point.
   *
   * @param time The time or tau coordinate (in appropriate units, e.g.,
   * GeV^-1).
   * @param x The space coordinate in the x-direction.
   * @param y The space coordinate in the y-direction.
   * @param z The space coordinate in the z-direction or eta coordinate.
   *
   * @return The temperature at the specified spacetime point in GeV.
   *
   * @see GetHydroInfo()
   */
  virtual Jetscape::real GetTemperature(Jetscape::real time, Jetscape::real x,
                                        Jetscape::real y, Jetscape::real z);

  /**
   * @brief Returns the fraction of quark-gluon plasma (QGP) at a given
   * spacetime point.
   *
   * This function computes the fraction of quark-gluon plasma assuming the
   * medium is in the QGP+Hadron Resonance Gas (HRG) phase at the specified
   * spacetime point (time, x, y, z).
   *
   * @param time The time or proper time (tau) coordinate.
   * @param x The space coordinate in the x-direction.
   * @param y The space coordinate in the y-direction.
   * @param z The space coordinate in the z-direction (or pseudorapidity eta).
   *
   * @return The fraction of quark-gluon plasma at the specified spacetime point
   * assuming medium is in QGP+HRG phase.
   *
   * @see GetHydroInfo()
   */
  virtual Jetscape::real GetQgpFraction(Jetscape::real time, Jetscape::real x,
                                        Jetscape::real y, Jetscape::real z);

  // These have no default implementation
  // /** @return 3-component (vx,vy,vz) fluid velocity at point (t or tau, x, y,
  // z or eta).
  //     @param time Time or tau coordinate.
  //     @param x Space coordinate.
  //     @param y Space coordinate.
  //     @param z Space or eta coordinate.
  // */
  // virtual real3 Get3FluidVelocity(Jetscape::real time, Jetscape::real x,
  // Jetscape::real y, Jetscape::real z)=0;

  // /** @return 4-component fluid velocity at point (t or tau, x, y, zor eta).
  //     @param time Time or tau coordinate.
  //     @param x Space coordinate.
  //     @param y Space coordinate.
  //     @param z Space or eta coordinate.
  // */
  // virtual real4 Get4FluidVelocity(Jetscape::real time, Jetscape::real x,
  // Jetscape::real y, Jetscape::real z);

  // /** @return Net baryon density at point (t or tau, x, y, z or eta).
  //     @param time Time or tau coordinate.
  //     @param x Space coordinate.
  //     @param y Space coordinate.
  //     @param z Space or eta coordinate.
  // */
  // virtual Jetscape::real GetNetBaryonDensity(Jetscape::real time,
  // Jetscape::real x, Jetscape::real y, Jetscape::real z);

  // /** @return Net charge density at point (t or tau, x, y, z or eta).
  // 	@param time Time or tau coordinate.
  // 	@param x Space coordinate.
  // 	@param y Space coordinate.
  // 	@param z Space or eta coordinate.
  // */
  // virtual Jetscape::real GetNetChargeDensity(Jetscape::real time,
  // Jetscape::real x, Jetscape::real y, Jetscape::real z);

  /**
   * @brief Adds a new liquefier to the system.
   *
   * This function accepts a shared pointer to a `LiquefierBase` object and
   * assigns it to the internal `liquefier_ptr`.
   *
   * @param new_liquefier A shared pointer to a `LiquefierBase` object to be
   * added.
   */
  virtual void add_a_liquefier(std::shared_ptr<LiquefierBase> new_liquefier) {
    liquefier_ptr = new_liquefier;
  }

  /**
   * @brief Retrieves the source term from the liquefier.
   *
   * This function locks the liquefier pointer and retrieves the source term
   * based on the given parameters. The source term is calculated and stored
   * in the provided array `jmu`.
   *
   * @param tau The proper time, typically used in relativistic fluid dynamics.
   * @param x The spatial coordinate in the x direction.
   * @param y The spatial coordinate in the y direction.
   * @param eta The spatial coordinate in the eta direction.
   * @param jmu An array of size 4 to store the resulting source term. The array
   *            will be populated with the computed values corresponding to the
   *            four components of the source term.
   *
   * @note The function relies on a locked pointer to the liquefier to access
   *       the source term calculation. This should be called when the liquefier
   *       object is initialized and the necessary resources are available.
   */
  void get_source_term(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                       Jetscape::real eta,
                       std::array<Jetscape::real, 4> jmu) const;

  /**
   * @brief Updates the energy deposit for a given time step.
   *
   * This function is intended to handle the "jet" signals in the future. It
   * logs the received jet signal with the time step `t` and energy deposit
   * value `edop` for debugging purposes.
   *
   * @param t The time step for the jet signal.
   * @param edop The energy deposit value associated with the jet signal.
   *
   * @note The function includes a commented-out line for a lock to ensure
   * thread safety (in multi-threaded environments).
   *
   * @todo Implement future functionality related to jet signals. TODO
   */
  virtual void UpdateEnergyDeposit(int t, double edop);

  /**
   * @brief Retrieves the energy density for the given time step.
   *
   * This function is a placeholder for future implementation of jet energy
   * density calculations. The energy density is set to zero by default.
   *
   * @param t The time step for which the energy density is to be retrieved.
   * @param edensity Reference to a variable where the energy density will be
   * stored.
   *
   * @note Currently, this function does not perform any calculations and sets
   * the energy density to 0.0.
   */
  virtual void GetEnergyDensity(int t, double &edensity) { edensity = 0.0; }

  /**
   * @brief Gets a reference to the bulk_info object.
   *
   * This function provides read-only access to the `bulk_info` object, which
   * contains data regarding the evolution history. It returns a constant
   * reference, ensuring that the caller cannot modify the object directly.
   *
   * @return A constant reference to the `EvolutionHistory` object.
   */
  const EvolutionHistory &get_bulk_info() const { return bulk_info; }

};  // end class FluidDynamics

}  // end namespace Jetscape

#endif  // FLUIDDYNAMICS_H
