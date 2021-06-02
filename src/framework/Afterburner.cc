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
// This is a general basic class for hadronic afterburner

#include "./Afterburner.h"
#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {
void Afterburner::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::Init();
  JSINFO << "Initializing Afterburner : " << GetId() << " ...";

  // Get the pointer to sampler
  soft_particlization_sampler_ =
      JetScapeSignalManager::Instance()->GetSoftParticlizationPointer().lock();
  if (!soft_particlization_sampler_) {
    JSWARN << "No soft particlization module found. It is necessary to provide"
           << " hadrons to afterburner.";
  }
  InitTask();
}

void Afterburner::Exec() {
  VERBOSE(2) << "Afterburner running: " << GetId() << " ...";
  ExecuteTask();
}

void Afterburner::CalculateTime() {
  VERBOSE(2) << "Afterburner running for time: " << GetId() << " ...";
  CalculateTimeTask();
}
} // end namespace Jetscape
