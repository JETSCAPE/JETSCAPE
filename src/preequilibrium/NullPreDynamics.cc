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

#include "NullPreDynamics.h"

#include <cstring>
#include <stdio.h>
#include <sys/stat.h>

#include "JetScapeLogger.h"

// Register the module with the base class
RegisterJetScapeModule<NullPreDynamics> NullPreDynamics::reg("NullPreDynamics");

NullPreDynamics::NullPreDynamics() {
  preequilibrium_status_ = NOT_STARTED;
  SetId("NullPreDynamics");
}

void NullPreDynamics::EvolvePreequilibrium() {
  VERBOSE(2) << "Initialize energy density profile in NullPreDynamics ...";
  // grab initial energy density from vector from initial state module
  std::vector<double> energy_density = ini->GetEntropyDensityDistribution();
  preequilibrium_status_ = INIT;
  if (preequilibrium_status_ == INIT) {
    VERBOSE(2) << "running NullPreDynamics ...";
    for (auto const &e_local : energy_density) {
      e_.push_back(e_local);
      P_.push_back(e_local / 3.);
      utau_.push_back(1.);
      ux_.push_back(0.);
      uy_.push_back(0.);
      ueta_.push_back(0.);
      pi00_.push_back(0.);
      pi01_.push_back(0.);
      pi02_.push_back(0.);
      pi03_.push_back(0.);
      pi11_.push_back(0.);
      pi12_.push_back(0.);
      pi13_.push_back(0.);
      pi22_.push_back(0.);
      pi23_.push_back(0.);
      pi33_.push_back(0.);
      bulk_Pi_.push_back(0.);
    }
    preequilibrium_status_ = DONE;
  }
}
