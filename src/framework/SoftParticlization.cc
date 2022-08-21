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

#include "SoftParticlization.h"

namespace Jetscape {

SoftParticlization::SoftParticlization() { boost_invariance = false; }

SoftParticlization::~SoftParticlization() {
  for (unsigned i = 0; i < Hadron_list_.size(); i++) {
    Hadron_list_.at(i).clear();
  }
  Hadron_list_.clear();
}

void SoftParticlization::Init() {
  JetScapeModuleBase::Init();
  JSINFO << "Initialize Soft particlization module ... " << GetId() << " ...";

  boost_invariance = check_boost_invariance();

  JSINFO << "boost invariance: " << boost_invariance;

  InitTask();
}

void SoftParticlization::Exec() {}

void SoftParticlization::Clear() {
  for (unsigned i = 0; i < Hadron_list_.size(); i++) {
    Hadron_list_.at(i).clear();
  }
  Hadron_list_.clear();
}

bool SoftParticlization::check_boost_invariance() {
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
