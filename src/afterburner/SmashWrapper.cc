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
// -----------------------------------------
// This is a wrapper for SMASH hadronic afterburner with the JETSCAPE framework
// -----------------------------------------

#include "SmashWrapper.h"

#include <filesystem>
#include <string>

#include "smash/library.h"
#include "smash/particles.h"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<SmashWrapper> SmashWrapper::reg("SMASH");

SmashWrapper::SmashWrapper() { SetId("SMASH"); }

void SmashWrapper::InitTask() {
  JSINFO << "SMASH: picking SMASH-specific configuration from xml file";
  std::string smash_config_file =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_config_file"});
  std::string smash_hadron_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_particles_file"});
  std::string smash_decays_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_decaymodes_file"});

  // do not store tabulation, which is achieved by an empty tabulations path
  std::string tabulations_path("");

  auto config = smash::setup_config_and_logging(
      smash_config_file, smash_hadron_list, smash_decays_list);

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config.set_value({"General", "Randomseed"}, random_seed);

  // Read in the rest of configuration
  end_time_ = GetXMLElementDouble({"Afterburner", "SMASH", "end_time"});
  config.set_value({"General", "End_Time"}, end_time_);

  JSINFO << "End time until which SMASH propagates is " << end_time_ << " fm/c";
  only_final_decays_ =
      GetXMLElementInt({"Afterburner", "SMASH", "only_decays"});
  if (only_final_decays_) {
    JSINFO << "SMASH will only perform resonance decays, no propagation";
  }

  const std::string smash_version(SMASH_VERSION);
  smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                     tabulations_path);

  JSINFO << "Seting up SMASH Experiment object";
  // output path is just dummy here, because no output from SMASH is foreseen
  std::filesystem::path output_path("./smash_output");
  smash_experiment_ =
      make_shared<smash::Experiment<AfterburnerModus>>(config, output_path);
  JSINFO << "Finish initializing SMASH";
}

void SmashWrapper::ExecuteTask() {
  AfterburnerModus *modus = smash_experiment_->modus();
  // This is necessary to correctly handle indices of particle sets from hydro.
  // Every hydro event creates a new structure like jetscape_hadrons_
  // with as many events in it as one has samples per hydro
  modus->reset_event_numbering();
  modus->jetscape_hadrons_ = GatherAfterburnerHadrons();
  const int n_events = modus->jetscape_hadrons_.size();
  JSINFO << "SMASH: obtained " << n_events << " events from particlization";
  // SMASH within JETSCAPE only works with one (the first) ensemble
  smash::Particles *smash_particles = smash_experiment_->first_ensemble();
  for (unsigned int i = 0; i < n_events; i++) {
    JSINFO << "Event " << i << " SMASH starts with "
           << modus->jetscape_hadrons_[i].size() << " particles.";
    smash_experiment_->initialize_new_event();
    if (!only_final_decays_) {
      smash_experiment_->run_time_evolution(end_time_);
    }
    smash_experiment_->do_final_decays();
    smash_experiment_->final_output();
    smash_particles_to_JS_hadrons(*smash_particles,
                                  modus->jetscape_hadrons_[i]);
    JSINFO << modus->jetscape_hadrons_[i].size() << " hadrons from SMASH.";
    smash_experiment_->increase_event_number();
  }
}

void SmashWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSINFO << "SMASH hadronic afterburner printout";
  auto f = w.lock();
  if (!f) {
    return;
  }
  AfterburnerModus *modus = smash_experiment_->modus();
  f->WriteComment("JetScape module: " + GetId());
  for (const auto &event : modus->jetscape_hadrons_) {
    int i = -1;
    for (const auto hadron : event) {
      f->WriteWhiteSpace("[" + to_string(++i) + "] H");
      f->Write(hadron);
    }
  }
}

void AfterburnerModus::JS_hadrons_to_smash_particles(
    const std::vector<shared_ptr<Hadron>> &JS_hadrons,
    smash::Particles &smash_particles) {
  smash_particles.reset();
  for (const auto JS_hadron : JS_hadrons) {
    const FourVector p = JS_hadron->p_in();
    const FourVector r = JS_hadron->x_in();
    const double mass = JS_hadron->restmass();
    smash::PdgCode pdgcode = smash::PdgCode::from_decimal(JS_hadron->pid());
    this->try_create_particle(smash_particles, pdgcode, r.t(), r.x(), r.y(),
                              r.z(), mass, JS_hadron->e(), JS_hadron->px(),
                              JS_hadron->py(), JS_hadron->pz());
  }
}

void SmashWrapper::smash_particles_to_JS_hadrons(
    const smash::Particles &smash_particles,
    std::vector<shared_ptr<Hadron>> &JS_hadrons) {
  JS_hadrons.clear();
  for (const auto &particle : smash_particles) {
    const int hadron_label = 0;
    const int hadron_status = 27;
    const int hadron_id = particle.pdgcode().get_decimal();
    smash::FourVector p = particle.momentum(), r = particle.position();
    const FourVector hadron_p(p.x1(), p.x2(), p.x3(), p.x0()),
        hadron_r(r.x1(), r.x2(), r.x3(), r.x0());
    const double hadron_mass = p.abs();
    JS_hadrons.push_back(make_shared<Hadron>(hadron_label, hadron_id,
                                             hadron_status, hadron_p, hadron_r,
                                             hadron_mass));
  }
}
