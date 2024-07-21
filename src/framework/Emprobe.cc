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
// -----------------------------------------
// JETSCPAE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// -----------------------------------------

#include "Emprobe.h"
#include "JetScapeSignalManager.h"

namespace Jetscape {

 Emprobe::Emprobe() {
 boost_invariance = false;
}


Emprobe::~Emprobe() {
}



void Emprobe::Init() {
  JetScapeModuleBase::Init();
  JSINFO << "Initialize EM probe module ... " << GetId() << " ...";

  boost_invariance = check_boost_invariance();

  JSINFO << "boost invariance: " << boost_invariance;

  InitTask();
}

void Emprobe::Exec() {}

void Emprobe::Clear() {
  
}

bool Emprobe::check_boost_invariance() {
  bool boost_invariance_flag = false;
  double grid_max_z = GetXMLElementDouble({"IS", "grid_max_z"});
  double grid_step_z = GetXMLElementDouble({"IS", "grid_step_z"});
  int nz = static_cast<int>(2. * grid_max_z / grid_step_z);
  if (nz <= 1) {
    boost_invariance_flag = true;
  }
  return (boost_invariance_flag);
}

} // end namespace Jetscape
