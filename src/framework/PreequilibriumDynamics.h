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

#ifndef PREEQUILDYNAMICS_H
#define PREEQUILDYNAMICS_H

#include <vector>

#include "FluidCellInfo.h"
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "RealType.h"

namespace Jetscape {

/// Enum for the status of the preequilibrium dynamics
enum PreequilibriumStatus { NOT_STARTED, INIT, DONE, ERR };

// Interface for the Preequilibrium Dynamics of the medium

/**
 * @brief JETSCAPE module for preequilibrium dynamics
 * 
 * This module will evolve the preequilibrium dynamics of the medium.
 */
class PreequilibriumDynamics : public JetScapeModuleBase {
 private:

 public:
  /**
   * @brief Construct a new PreequilibriumDynamics object
   */
  PreequilibriumDynamics();

  /**
   * @brief Destroy the PreequilibriumDynamics object
   */
  virtual ~PreequilibriumDynamics();

  /// Initial and final time for preequilibrium evolution
  real preequilibrium_tau_0_, preequilibrium_tau_max_;

  /**
   * @brief Initialize the PreequilibriumDynamics module
   * 
   * Reads the input parameters from the XML file under the tag <Preequilibrium>.
   * Uses JetScapeSingnalManager instance to retrieve the Initial State Physics
   * information.
   */
  void Init();
  
  /**
   * @brief Execute the PreequilibriumDynamics module
   * 
   * Calls EvolvePreequilibrium(). This explicit call can be used for actual
   * execution of pre-equilibrium evolution defined in the modules.
   */
  void Exec();

  /**
   * @brief Clear the PreequilibriumDynamics module
   * 
   * This function clears the preequilibrium dynamics module. It can be
   * overridden by other tasks.
   */
  virtual void Clear();

  /**
   * @brief Initialize the preequilibrium dynamics
   * 
   * This function initializes the preequilibrium dynamics. It can be overridden
   * by other tasks.
   */
  virtual void InitializePreequilibrium() {}
  
  /**
   * @brief Evolve the preequilibrium dynamics
   * 
   * This function evolves the preequilibrium dynamics. It can be overridden by
   * other tasks.
   */
  virtual void EvolvePreequilibrium() {}

  /**
   * @brief Set the initial state object
   * 
   * @param ini A shared pointer of type InitialState
   */
  std::shared_ptr<InitialState> ini;

  /**
   * @brief Get the pre-equilibrium status
   * 
   * @return int The status of the preequilibrium dynamics 
   */
  int GetPreequilibriumStatus() { return (preequilibrium_status_); }

  // @return Start time (or tau) for hydrodynamic evolution

  /**
   * @brief Get the preequilibrium start time
   * 
   * @return real The start time for preequilibrium evolution
   */
  virtual real GetPreequilibriumStartTime() const {
    return (preequilibrium_tau_0_);
  }

  /**
   * @brief Get the preequilibrium evolution time step
   * 
   * This function can be overridden by other tasks.
   * 
   * @return real The time step for preequilibrium evolution
   */
  virtual real GetPreequilibriumEvodtau() const { return (0.02); }

  /**
   * @brief Get the number of time steps for preequilibrium evolution
   * 
   * This function can be overridden by other tasks.
   * 
   * @return int The number of time steps for preequilibrium evolution
   */
  virtual int get_ntau() const { return (0); }

  /**
   * @brief Get the preequilibrium end time
   * 
   * @return real The end time for preequilibrium evolution
   */
  real GetPreequilibriumEndTime() { return (preequilibrium_tau_max_); }

  /**
   * @brief Get the number of fluid cells
   * 
   * This function can be overridden by other tasks.
   * 
   * @return int The number of fluid cells
   */
  virtual int get_number_of_fluid_cells() { return (0); }
  
  /**
   * @brief Get the fluid cell with index
   * 
   * This function can be overridden by other tasks.
   * 
   * @param idx The index of the fluid cell
   * @param info_ptr A unique pointer of type FluidCellInfo
   */
  virtual void get_fluid_cell_with_index(
      const int idx, std::unique_ptr<FluidCellInfo> &info_ptr) {}

  /**
   * @brief Clear the evolution data
   * 
   * This function can be overridden by other tasks.
   */
  virtual void clear_evolution_data() {}

  /// Pre-equilibrium status
  PreequilibriumStatus preequilibrium_status_;

  /// Pre-equilibrium vectors containing the bulk properties
  std::vector<double> e_;
  std::vector<double> P_;
  std::vector<double> utau_;
  std::vector<double> ux_;
  std::vector<double> uy_;
  std::vector<double> ueta_;
  std::vector<double> pi00_;
  std::vector<double> pi01_;
  std::vector<double> pi02_;
  std::vector<double> pi03_;
  std::vector<double> pi11_;
  std::vector<double> pi12_;
  std::vector<double> pi13_;
  std::vector<double> pi22_;
  std::vector<double> pi23_;
  std::vector<double> pi33_;
  std::vector<double> bulk_Pi_;
};

}  // end namespace Jetscape

#endif
