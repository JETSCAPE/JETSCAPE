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
// This is a general basic class for hadronic afterburner

#include "./Afterburner.h"

#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {
void Afterburner::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::Init();
  JSINFO << "Initializing Afterburner : " << GetId() << " ...";
  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};
  InitTask();
}

void Afterburner::Exec() {
  VERBOSE(2) << "Afterburner running: " << GetId() << " ...";
  ExecuteTask();
}

std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GetSoftParticlizationHadrons() {
  auto soft_particlization =
      JetScapeSignalManager::Instance()->GetSoftParticlizationPointer().lock();
  if (!soft_particlization) {
    JSWARN << "No soft particlization module found. Check if fragmentation"
           << " hadrons are handed to afterburner.";
    std::vector<std::shared_ptr<Hadron>> hadrons;
    dummy.push_back(hadrons);
    return dummy;
  } else {
    return soft_particlization->Hadron_list_;
  }
}

std::vector<shared_ptr<Hadron>> Afterburner::GetFragmentationHadrons() {
  JSINFO << "Get fragmentation hadrons in Afterburner";
  auto hadronization_mgr = JetScapeSignalManager::Instance()
                               ->GetHadronizationManagerPointer()
                               .lock();
  if (!hadronization_mgr) {
    JSWARN << "No hardronization module found. It is necessary to include"
           << " fragmentation hadrons to afterburner as requested.";
    exit(1);
  }
  std::vector<shared_ptr<Hadron>> h_list;
  hadronization_mgr->GetHadrons(h_list);
  JSINFO << "Got " << h_list.size()
         << " fragmentation hadrons from HadronizationManager.";

  std::vector<shared_ptr<Hadron>> h_list_new;
  rand_int_ptr_ = (std::make_shared<std::uniform_int_distribution<int>>(0, 1));
  for (auto h : h_list) {
    if (h->has_no_position()) {
      JSDEBUG << "Found fragmentation hadron without properly set position in "
                 "Afterburner.\nInclusion of fragmentation hadrons only "
                 "possible for HybridHadronization.";
    }

    // move all the fragmentation hadrons a little bit around to avoid having
    // multiple hadrons at the same position if they are at the same position
    const FourVector r = h->x_in();
    const double rand_x =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    const double rand_y =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    const double rand_z =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    double position_smeared[4] = {r.t(), r.x() + rand_x, r.y() + rand_y,
                                  r.z() + rand_z};
    h->set_x(position_smeared);

    if ((std::abs(h->pid()) > 10) && (h->pid() != 21)) {
      if (h->pstat() > 0) {
        // convert Kaon-L or Kaon-S into K0 or Anti-K0
        if (h->pid() == 310 || h->pid() == 130) {
          const int rand_int = (*rand_int_ptr_)(*GetMt19937Generator());
          const int id = (rand_int == 0) ? 311 : -311;
          h->set_id(id);
        }
        h_list_new.push_back(h);
      } else if (h->pstat() < 0) {
        // convert Kaon-L or Kaon-S into K0 or Anti-K0
        // change id of negative Kaons to make them consistent with the SMASH
        // output
        if (h->pid() == 310 || h->pid() == 130) {
          const int rand_int = (*rand_int_ptr_)(*GetMt19937Generator());
          const int id = (rand_int == 0) ? 311 : -311;
          h->set_id(id);
        }
      }
    } else if ((std::abs(h->pid()) < 10) || (h->pid() == 21)) {
      JSWARN << "Found a free quark or gluon! This can not be handed over to "
                "SMASH.\n"
                "Check confinement in hadronization module!";
    }
  }
  return h_list_new;
}

std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GatherAfterburnerHadrons() {
  std::vector<std::vector<shared_ptr<Hadron>>> afterburner_had_events;
  afterburner_had_events = GetSoftParticlizationHadrons();

  if (GetXMLElementInt({"Afterburner", "output_only_final_state_hadrons"})) {
    // clear Hadron_list_ in soft_particlization, otherwise the final hadron
    // output of the writer contains also the soft hadrons which were used as
    // input for SMASH
    auto soft_particlization = JetScapeSignalManager::Instance()
                                   ->GetSoftParticlizationPointer()
                                   .lock();
    if (soft_particlization) {
      soft_particlization->Hadron_list_.clear();
    }
  }

  if (GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"})) {
    if (afterburner_had_events.size() > 1) {
      JSWARN
          << "Fragmentation hadrons in Afterburner are only possible without "
             "repeated sampling from SoftParticlization. Exiting.";
      exit(1);
    }
    std::vector<shared_ptr<Hadron>> frag_hadrons = GetFragmentationHadrons();

    if (GetXMLElementInt({"Afterburner", "output_only_final_state_hadrons"})) {
      // empty the hadron vector in the hadronization manager to circumvent the
      // output of these hadrons if they are implemented in the SMASH
      // afterburner
      auto hadronization_mgr = JetScapeSignalManager::Instance()
                                   ->GetHadronizationManagerPointer()
                                   .lock();
      hadronization_mgr->DeleteRealHadrons();
    }

    afterburner_had_events[0].insert(afterburner_had_events[0].end(),
                                     frag_hadrons.begin(), frag_hadrons.end());
    dummy.clear();
  }
  return afterburner_had_events;
}

}  // end namespace Jetscape
