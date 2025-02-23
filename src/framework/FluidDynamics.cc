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

#include "FluidDynamics.h"

#include <array>
#include <iostream>

#include "JetScapeSignalManager.h"
#include "LinearInterpolation.h"
#include "MakeUniqueHelper.h"
#include "SurfaceFinder.h"

#define MAGENTA "\033[35m"

using namespace std;

namespace Jetscape {

/**
 * @brief Default constructor for the FluidDynamics class.
 *
 * Initializes the object with a default task ID of "FluidDynamics" and
 * sets the eta parameter to -99.99. Additionally, the boost_invariant_
 * is set to true. This constructor is typically used to initialize
 * an instance of FluidDynamics with default values.
 *
 * @note The task ID is set via the SetId method to "FluidDynamics".
 * @note The eta parameter is initialized to -99.99 to indicate an uninitialized
 * state.
 */
FluidDynamics::FluidDynamics() {
  VERBOSE(8);
  eta = -99.99;
  boost_invariant_ = true;
  SetId("FluidDynamics");
}

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
FluidDynamics::~FluidDynamics() {
  VERBOSE(8);
  disconnect_all();
}

/**
 * @brief Initializes the FluidDynamics module by reading input parameters
 *        from the XML file under the <Hydro> tag, retrieving the initial state
 *        physics information, and initializing related modules and tasks.
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
void FluidDynamics::Init() {
  JetScapeModuleBase::Init();

  JSINFO << "Initialize FluidDynamics : " << GetId() << " ...";

  VERBOSE(8);
  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
    JSWARN << "No initialization module, "
           << "try: auto trento = make_shared<TrentoInitial>(); "
           << "jetscape->Add(trento);";
  }

  pre_eq_ptr =
      JetScapeSignalManager::Instance()->GetPreEquilibriumPointer().lock();
  if (!pre_eq_ptr) {
    JSWARN << "No Pre-equilibrium module";
  }

  InitializeHydro(parameter_list);
  InitTask();

  JetScapeTask::InitTasks();
}

/**
 * @brief Executes the hydrodynamic evolution and additional tasks within the
 * current module.
 *
 * This function explicitly calls the `EvolveHydro()` method to perform the
 * hydrodynamic evolution, which is defined in various modules like @a Brick, @a
 * MpiMusic, or @a OSU-HYDRO if they are attached as polymorphic classes. In
 * addition, the function also executes the tasks associated with the current
 * module.
 *
 * It provides verbose logging for debugging purposes, including information
 * about the system state such as entropy density distribution size and event
 * number.
 *
 * @sa Read about @a polymorphism in C++ for more information on how polymorphic
 * classes are used in the context of hydrodynamic evolution.
 *
 * @note This function is designed to be virtual, allowing for potential
 * overrides in derived classes to implement specific behavior.
 *
 * @see EvolveHydro()
 * @see JetScapeTask::ExecuteTasks()
 */
void FluidDynamics::Exec() {
  VERBOSE(2) << "Run Hydro : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  if (ini) {
    VERBOSE(3) << "length of entropy density vector="
               << ini->GetEntropyDensityDistribution().size();
  }

  EvolveHydro();
  JetScapeTask::ExecuteTasks();
}

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
void FluidDynamics::Clear() {
  clear_up_evolution_data();
  if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
    liquefier_ptr.lock()->Clear();
  }
}

/**
 * @brief Collects header information for writer modules.
 *
 * This function extracts and sets the event plane angle in the header of the
 * provided JetScapeWriter instance, if the writer is valid.
 *
 * @param w A weak pointer to a JetScapeWriter instance. This pointer is used to
 *          retrieve the header for setting the event plane angle.
 *
 * @note If the weak pointer is expired or invalid, no operation is performed.
 */
void FluidDynamics::CollectHeader(weak_ptr<JetScapeWriter> w) {
  auto f = w.lock();
  if (f) {
    auto &header = f->GetHeader();
    header.SetEventPlaneAngle(GetEventPlaneAngle());
  }
}

/**
 * @brief Calculates the hypersurface for the Cooper-Frye or recombination
 * model.
 *
 * This function is designed to return the hypersurface at a given constant
 * temperature and can be overridden by different modules for more specific
 * implementations. It is part of the fluid dynamics computations, where the
 * hypersurface is needed to track the evolution of the system at the specified
 * temperature.
 *
 * @param T_sw The switching temperature at which the hypersurface is computed.
 *             This is a scalar value representing the temperature threshold.
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
void FluidDynamics::FindAConstantTemperatureSurface(
    Jetscape::real T_sw, std::vector<SurfaceCellInfo> &surface_cells) {
  std::unique_ptr<SurfaceFinder> surface_finder_ptr(
      new SurfaceFinder(T_sw, bulk_info));
  surface_finder_ptr->Find_full_hypersurface();
  surface_cells = surface_finder_ptr->get_surface_cells_vector();
  JSINFO << "number of surface cells: " << surface_cells.size();
}

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
Jetscape::real FluidDynamics::GetEnergyDensity(Jetscape::real time,
                                               Jetscape::real x,
                                               Jetscape::real y,
                                               Jetscape::real z) {
  std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
  GetHydroInfo(time, x, y, z, fluid_cell_ptr);
  real energy_density = fluid_cell_ptr->energy_density;
  return (energy_density);
}

/**
 * @brief Returns the entropy density at a given space-time point.
 *
 * This function computes and returns the entropy density at a specified
 * space-time point (time, x, y, z). It retrieves the fluid cell information
 * from the given space-time coordinates and accesses the entropy density value
 * stored within the fluid cell.
 *
 * @param time The time or tau coordinate of the space-time point.
 * @param x The x-space coordinate of the space-time point.
 * @param y The y-space coordinate of the space-time point.
 * @param z The z-space coordinate or eta coordinate of the space-time point.
 *
 * @return The entropy density at the given space-time point, in units of [GeV].
 *
 * @see GetHydroInfo()
 */
Jetscape::real FluidDynamics::GetEntropyDensity(Jetscape::real time,
                                                Jetscape::real x,
                                                Jetscape::real y,
                                                Jetscape::real z) {
  std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
  GetHydroInfo(time, x, y, z, fluid_cell_ptr);
  real entropy_density = fluid_cell_ptr->entropy_density;
  return (entropy_density);
}

/**
 * @brief Returns the temperature at a specific spacetime point.
 *
 * This function calculates and returns the temperature at the given spacetime
 * coordinates (time, x, y, z) in GeV. It retrieves the necessary fluid dynamics
 * information and extracts the temperature at the specified point.
 *
 * @param time The time or tau coordinate (in appropriate units, e.g., GeV^-1).
 * @param x The space coordinate in the x-direction.
 * @param y The space coordinate in the y-direction.
 * @param z The space coordinate in the z-direction or eta coordinate.
 *
 * @return The temperature at the specified spacetime point in GeV.
 *
 * @see GetHydroInfo()
 */
Jetscape::real FluidDynamics::GetTemperature(Jetscape::real time,
                                             Jetscape::real x, Jetscape::real y,
                                             Jetscape::real z) {
  std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
  GetHydroInfo(time, x, y, z, fluid_cell_ptr);
  real temperature = fluid_cell_ptr->temperature;
  return (temperature);
}

/**
 * @brief Returns the fraction of quark-gluon plasma (QGP) at a given spacetime
 * point.
 *
 * This function computes the fraction of quark-gluon plasma assuming the medium
 * is in the QGP+Hadron Resonance Gas (HRG) phase at the specified spacetime
 * point (time, x, y, z).
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
Jetscape::real FluidDynamics::GetQgpFraction(Jetscape::real time,
                                             Jetscape::real x, Jetscape::real y,
                                             Jetscape::real z) {
  std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
  GetHydroInfo(time, x, y, z, fluid_cell_ptr);
  real qgp_fraction = fluid_cell_ptr->qgp_fraction;
  return (qgp_fraction);
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
void FluidDynamics::get_source_term(Jetscape::real tau, Jetscape::real x,
                                    Jetscape::real y, Jetscape::real eta,
                                    std::array<Jetscape::real, 4> jmu) const {
  liquefier_ptr.lock()->get_source(tau, x, y, eta, jmu);
}

/**
 * @brief Prints out the information of the fluid cell to the screen.
 *
 * This function outputs detailed information about the fluid cell's state,
 * including various physical properties such as energy density, entropy
 * density, temperature, pressure, chemical potentials, velocity components, and
 * shear/bulk viscosities.
 *
 * @param fluid_cell_info_ptr A pointer to the `FluidCellInfo` object containing
 * the fluid cell's data.
 *
 * @note This function outputs the information using the `JSINFO` logging
 * system. The values printed are in appropriate units such as GeV/fm^3 for
 * energy and pressure, and GeV for chemical potentials.
 *
 * @see FluidCellInfo
 */
void FluidDynamics::PrintFluidCellInformation(
    FluidCellInfo *fluid_cell_info_ptr) {
  // this function print out the information of the fluid cell to the screen
  JSINFO << "=======================================================";
  JSINFO << "print out cell information:";
  JSINFO << "=======================================================";
  JSINFO << "energy density = " << fluid_cell_info_ptr->energy_density
         << " GeV/fm^3.";
  JSINFO << "entropy density = " << fluid_cell_info_ptr->entropy_density
         << " 1/fm^3.";
  JSINFO << "temperature = " << fluid_cell_info_ptr->temperature << " GeV.";
  JSINFO << "pressure = " << fluid_cell_info_ptr->pressure << " GeV/fm^3.";
  JSINFO << "QGP_fraction = " << fluid_cell_info_ptr->qgp_fraction;
  JSINFO << "mu_B = " << fluid_cell_info_ptr->mu_B << " GeV.";
  JSINFO << "mu_S = " << fluid_cell_info_ptr->mu_S << " GeV.";
  JSINFO << "mu_C = " << fluid_cell_info_ptr->mu_C << " GeV.";
  JSINFO << "vx = " << fluid_cell_info_ptr->vx;
  JSINFO << "vy = " << fluid_cell_info_ptr->vy;
  JSINFO << "vz = " << fluid_cell_info_ptr->vz;
  JSINFO << "shear viscous pi^{munu} (GeV/fm^3): ";
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      JSINFO << fluid_cell_info_ptr->pi[i][j];
    }
  }
  JSINFO << "bulk_Pi = " << fluid_cell_info_ptr->bulk_Pi << " GeV/fm^3";
  JSINFO << "=======================================================";
}

/**
 * @brief Updates the energy deposit for a given time step.
 *
 * This function is intended to handle the "jet" signals in the future. It logs
 * the received jet signal with the time step `t` and energy deposit value
 * `edop` for debugging purposes.
 *
 * @param t The time step for the jet signal.
 * @param edop The energy deposit value associated with the jet signal.
 *
 * @note The function includes a commented-out line for a lock to ensure thread
 * safety (in multi-threaded environments).
 *
 * @todo Implement future functionality related to jet signals. TODO
 */
void FluidDynamics::UpdateEnergyDeposit(int t, double edop) {
  // sigslot::lock_block<multi_threaded_local> lock(this);
  JSDEBUG << MAGENTA << "Jet Signal received : " << t << " " << edop;
}

}  // end namespace Jetscape
