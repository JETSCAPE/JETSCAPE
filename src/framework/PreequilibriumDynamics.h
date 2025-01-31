/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "FluidCellInfo.h"
#include "RealType.h"

namespace Jetscape {
// Flags for preequilibrium dynamics status.
enum PreequilibriumStatus { NOT_STARTED, INIT, DONE, ERR };


// Interface for the Preequilibrium Dynamics of the medium
class PreequilibriumDynamics : public JetScapeModuleBase {
private:

public:
  PreequilibriumDynamics();

  virtual ~PreequilibriumDynamics();
    real preequilibrium_tau_0_, preequilibrium_tau_max_;

  /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls InitializeHydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
  void Init();

  /** Calls EvolvePreequilibrium(); This explicit call can be used for actual execution of Preequilibrium evolution defined in the modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
  void Exec();

  /** Default Clear() function. It can be overridden by other tasks.*/
  virtual void Clear();

  virtual void InitializePreequilibrium() {}
  virtual void EvolvePreequilibrium() {}

  // add initial state shared pointer
  /** A pointer of type InitialState class.
    */
  std::shared_ptr<InitialState> ini;

  int GetPreequilibriumStatus() { return (preequilibrium_status_); }

  // @return Start time (or tau) for hydrodynamic evolution
  virtual real GetPreequilibriumStartTime() const {
      return (preequilibrium_tau_0_);
  }

  virtual real GetPreequilibriumEvodtau() const { return (0.02); }

  virtual int get_ntau() const { return(0); }

  // @return End time (or tau) for hydrodynamic evolution.
  real GetPreequilibriumEndTime() { return (preequilibrium_tau_max_); }

  virtual int get_number_of_fluid_cells() { return(0); }
  virtual void get_fluid_cell_with_index(
          const int idx, std::unique_ptr<FluidCellInfo> &info_ptr) {}
  virtual void clear_evolution_data() {}

  // record preequilibrium running status
  PreequilibriumStatus preequilibrium_status_;

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

} // end namespace Jetscape

#endif
