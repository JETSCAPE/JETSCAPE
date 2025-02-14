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

#ifndef FREESTREAMMILNEWRAPPER_H
#define FREESTREAMMILNEWRAPPER_H

#include "FreestreamMilne.cpp"
#include "PreequilibriumDynamics.h"

using namespace Jetscape;

class FreestreamMilneWrapper : public PreequilibriumDynamics {
 private:
  // int mode; //!< records running mode
  FREESTREAMMILNE *fsmilne_ptr;

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<FreestreamMilneWrapper> reg;

 public:
  FreestreamMilneWrapper();
  ~FreestreamMilneWrapper();

  void InitializePreequilibrium();
  void EvolvePreequilibrium();

  int get_number_of_fluid_cells() {
    return (fsmilne_ptr->get_number_of_fluid_cells());
  }

  void get_fluid_cell_with_index(const int idx,
                                 std::unique_ptr<FluidCellInfo> &info_ptr);
  void clear_evolution_data() { fsmilne_ptr->clear_evolution_data(); }

  Jetscape::real GetPreequilibriumStartTime() const {
    return (fsmilne_ptr->GetPreequilibriumStartTime());
  }

  Jetscape::real GetPreequilibriumEvodtau() const {
    return (fsmilne_ptr->GetPreequilibriumEvodtau());
  }

  int get_ntau() const { return (fsmilne_ptr->get_ntau()); }
};

#endif  // FREESTREAMMILNEWRAPPER_H
