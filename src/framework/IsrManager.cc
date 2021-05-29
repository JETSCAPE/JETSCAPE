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

#include "IsrManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "MakeUniqueHelper.h"
#include <string>

#include <iostream>
#include <vector>
#include <thread>

using namespace std;

namespace Jetscape {

IsrManager::IsrManager() : JetEnergyLossManager() {
  SetId("IsrManager");
  VERBOSE(8);
}

IsrManager::~IsrManager() {
  // Check if this is all really needed with shared_ptr ...
  JSDEBUG;
  Clear();

  if (GetNumberOfTasks() > 0)
    EraseTaskLast();
}

void IsrManager::Init()
{
  JSINFO << "Intialize ISR Manager ..."; //via JetEnergyLossManager::Init() ...";

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid ISR Manager modules found ...";
    exit(-1);
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " ISR Manager Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  JSINFO << "Connect IsrManager Signal to Hard Process ...";
  JetScapeSignalManager::Instance()->ConnectGetHardPartonListSignal(
      dynamic_pointer_cast<IsrManager>(shared_from_this()));

  /*
  // Set the pointer of JetEnergyLoss for making connections to hadronization module
  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<JetEnergyLoss>(it))

      JetScapeSignalManager::Instance()->SetEnergyLossPointer(
          dynamic_pointer_cast<JetEnergyLoss>(it));
  }
  */

  //JetEnergyLossManager::Init();
}

void IsrManager::Exec()
{
  VERBOSE(1) << "Run ISR Manager ...";
  JSDEBUG << "Task Id = " << this_thread::get_id();

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid ISR Manager modules found ...";
    exit(-1);
  }

  JetEnergyLossManager::Exec();
}

} // end namespace Jetscape
