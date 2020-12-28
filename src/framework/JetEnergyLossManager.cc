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

#include "JetEnergyLossManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "MakeUniqueHelper.h"
#include <string>

#include <iostream>
#include <vector>
#include <thread>

using namespace std;

namespace Jetscape {

JetEnergyLossManager::JetEnergyLossManager() {
  SetId("JLossManager");
  GetHardPartonListConnected = false;
  VERBOSE(8);
}

JetEnergyLossManager::~JetEnergyLossManager() {
  // Check if this is all really needed with shared_ptr ...
  JSDEBUG;
  Clear();

  if (GetNumberOfTasks() > 0)
    EraseTaskLast();
}

void JetEnergyLossManager::Clear() {
  JSDEBUG << "Hard Parton List ...";

  hp.clear();

  int n = GetNumberOfTasks();
  for (int i = 1; i < n; i++)
    EraseTaskLast();

  // Clean Up not really working with iterators (see also above!!!) Some logic not clear for me.
  JetScapeSignalManager::Instance()->CleanUp();
  JetScapeTask::ClearTasks();

  VERBOSE(8) << hp.size();
}

void JetEnergyLossManager::Init() {
  JSINFO << "Intialize JetEnergyLoss Manager ...";

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Energy Loss Manager modules found ...";
    exit(-1);
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " Eloss Manager Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  JSINFO << "Connect JetEnergyLossManager Signal to Hard Process ...";
  JetScapeSignalManager::Instance()->ConnectGetHardPartonListSignal(
      shared_from_this());

  // Set the pointer of JetEnergyLoss for making connections to hadronization module
  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<JetEnergyLoss>(it))

      JetScapeSignalManager::Instance()->SetEnergyLossPointer(
          dynamic_pointer_cast<JetEnergyLoss>(it));
  }
}

void JetEnergyLossManager::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLossManager::Exec() {
  VERBOSE(1) << "Run JetEnergyLoss Manager ...";
  JSDEBUG << "Task Id = " << this_thread::get_id();

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Energy Loss Manager modules found ...";
    exit(-1);
  }

  // ----------------------------------
  // Create needed copies and connect signal/slots accordingly ...

  if (GetGetHardPartonListConnected()) {
    GetHardPartonList(hp);
    VERBOSE(3) << " Number of Hard Partons = " << hp.size();
    for (int i = 1; i < hp.size(); i++) {
      JSDEBUG << "Create the " << i
              << " th copy because number of intital hard partons = "
              << hp.size();
      // Add(make_shared<JetEnergyLoss>(*dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0))));
      auto jloss_org = dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0));
      auto jloss_copy = make_shared<JetEnergyLoss>(*jloss_org);

      // if there is a liquefier attached to the jloss module
      // also attach the liquefier to the copied jloss modules
      // to collect hydrodynamic source terms
      if (!weak_ptr_is_uninitialized(jloss_org->get_liquefier())) {
        jloss_copy->add_a_liquefier(jloss_org->get_liquefier().lock());
      }
      Add(jloss_copy);
    }
  }

  VERBOSE(3) << " Found " << GetNumberOfTasks()
             << " Eloss Manager Tasks/Modules Execute them ... ";
  JSDEBUG << "Check and Create Signal/Slots via JetScapeSignalManaher instance "
             "if needed ...";

  CreateSignalSlots();

  // ----------------------------------
  // Copy Shower initiating partons to JetEnergyLoss tasks ...

  if (GetGetHardPartonListConnected() && hp.size() > 0) {
    int n = 0;
    for (auto it : GetTaskList()) {
      dynamic_pointer_cast<JetEnergyLoss>(it)->AddShowerInitiatingParton(
          hp.at(n));
      n++;
    }
  }

  // ----------------------------------
  // quick and dirty here, only include after further testing (flag in init xml files ...)
  bool multiTask = false;

  // ----------------------------------
  //Excute JetEnergyLoss tasks and their subtasks (done via signal/slot) by hand ...
  //needed if only JetEnergyloss tasks in parallel (otherwise could be done via JetScapeTask, then every task a new thread for example)
  // Of course now the number of threads if plainly given by the number of initial hard partons and can exceed the allowed # of threads >1000
  // so to really do this properly one has to think about and limit to number of CPU's * N threads or something in that directions. See a quick attempt below.

  //DEBUG:
  //unsigned num_cpus = thread::hardware_concurrency();
  //cout << "Num of CPU's = " << num_cpus << " threads\n";

  if (multiTask) {
    int nTasks = GetNumberOfTasks();
    int nCPUs = thread::hardware_concurrency();

    vector<thread> threads;

    int nMaxThreads = nCPUs * 2;
    int n = 0;

    VERBOSE(3) << " Use multi-threading: (max) # of threads = # of CPU's "
               << nCPUs << " (found) * 2";

    for (auto it : GetTaskList()) {
      if (it->GetActive()) {
        threads.push_back(thread(&JetEnergyLoss::Exec,
                                 dynamic_pointer_cast<JetEnergyLoss>(it)));
        n++;
      }
      if (n == nMaxThreads) {
        //DEBUGTHREAD<<n;
        for (auto &th : threads)
          th.join();
        n = 0;
        threads.clear();

        //DEBUG:
        //std::this_thread::sleep_for(std::chrono::milliseconds(5000));
      }
    }

    if (nTasks < nMaxThreads) {
      for (auto &th : threads)
        th.join();
      threads.clear();
    }
  }
  // ----------------------------------
  else
    // Standard "serial" execution for the JetEnerguLoss (+submodules) task ...
    JetScapeTask::ExecuteTasks();

  //Add acheck if the parton shower was actually created for the Modules ....
  VERBOSE(3) << " " << GetNumberOfTasks()
             << " Eloss Manager Tasks/Modules finished.";
}

void JetEnergyLossManager::CreateSignalSlots() {
  for (auto it : GetTaskList()) {
    for (auto it2 : it->GetTaskList()) {
      if (!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetJetSignalConnected()) {
        JetScapeSignalManager::Instance()->ConnectJetSignal(
            dynamic_pointer_cast<JetEnergyLoss>(it2));
      }
      if (!dynamic_pointer_cast<JetEnergyLoss>(it2)
               ->GetEdensitySignalConnected()) {
        JetScapeSignalManager::Instance()->ConnectEdensitySignal(
            dynamic_pointer_cast<JetEnergyLoss>(it2));
      }
      if (!dynamic_pointer_cast<JetEnergyLoss>(it2)
               ->GetGetHydroCellSignalConnected()) {
        JetScapeSignalManager::Instance()->ConnectGetHydroCellSignal(
            dynamic_pointer_cast<JetEnergyLoss>(it2));
      }
      if (!dynamic_pointer_cast<JetEnergyLoss>(it2)
               ->GetGetHydroTau0SignalConnected()) {
        JetScapeSignalManager::Instance()->ConnectGetHydroTau0Signal(
            dynamic_pointer_cast<JetEnergyLoss>(it2));
      }

      // between eloss modules and eloss
      // check the signals itself, probably best via manager in the long run ...
      if (!dynamic_pointer_cast<JetEnergyLoss>(it2)
               ->GetSentInPartonsConnected()) {
        JetScapeSignalManager::Instance()->ConnectSentInPartonsSignal(
            dynamic_pointer_cast<JetEnergyLoss>(it),
            dynamic_pointer_cast<JetEnergyLoss>(it2));
      }
    }
    auto liq_pt = dynamic_pointer_cast<JetEnergyLoss>(it)->get_liquefier();
    if (!weak_ptr_is_uninitialized(liq_pt) &&
        !liq_pt.lock()->get_GetHydroCellSignalConnected()) {
      JetScapeSignalManager::Instance()->ConnectGetHydroCellSignal(
          liq_pt.lock());
    }
  }

  JetScapeSignalManager::Instance()->PrintGetHydroCellSignalMap();
  VERBOSE(8);
  JetScapeSignalManager::Instance()->PrintSentInPartonsSignalMap();
}

} // end namespace Jetscape
