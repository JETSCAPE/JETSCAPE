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
// This is a wrapper for SMASH hadronic afterburner with the JETSCAPE framework
// -----------------------------------------

#include "SmashWrapper.h"

#include "smash/decaymodes.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"
#include "smash/sha256.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<SmashWrapper> SmashWrapper::reg("SMASH");

SmashWrapper::SmashWrapper() { SetId("SMASH"); }

void SmashWrapper::InitTask() {
  JSINFO << "SMASH: picking SMASH-specific configuration from xml file";
  std::string smash_config =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_config_file"});
  boost::filesystem::path input_config_path(smash_config);
  if (!boost::filesystem::exists(input_config_path)) {
    JSWARN << "SMASH config file " << smash_config << " not found.";
    exit(-1);
  } else {
    JSINFO << "Obtaining SMASH configuration from " << smash_config;
  }
  smash::Configuration config(input_config_path.parent_path(),
                              input_config_path.filename());
  // SMASH hadron list
  std::string hadron_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_particles_file"});
  if (boost::filesystem::exists(hadron_list)) {
    config["particles"] =
        smash::read_all(boost::filesystem::ifstream{hadron_list});
    smash::ParticleType::create_type_list(config.take({"particles"}));
    JSINFO << "Obtaining SMASH hadron list from " << hadron_list;
  } else {
    JSINFO << "Using default SMASH hadron list.";
  }
  // SMASH decaymodes
  std::string decays_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_decaymodes_file"});
  if (boost::filesystem::exists(decays_list)) {
    config["decaymodes"] =
        smash::read_all(boost::filesystem::ifstream{decays_list});
    smash::DecayModes::load_decaymodes(config.take({"decaymodes"}));
    JSINFO << "Obtaining SMASH decays list from " << decays_list;
  } else {
    JSINFO << "Using default SMASH decay list.";
  }
  // Make sure that unstable hadrons have decays and stable not
  smash::ParticleType::check_consistency();

  // Let SMASH tabulate cross sections
  boost::filesystem::path tabulations_path = "./tabulations";
  boost::filesystem::create_directories(tabulations_path);
  const std::string particle_string = config["particles"].to_string();
  const std::string decay_string = config["decaymodes"].to_string();
  smash::sha256::Context hash_context;
  hash_context.update(particle_string);
  hash_context.update(decay_string);
  const auto hash = hash_context.finalize();
  smash::IsoParticleType::tabulate_integrals(hash, tabulations_path);

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config["General"]["Randomseed"] = random_seed;
  // Set SMASH logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::TRACE));
  smash::create_all_loggers(config["Logging"]);
  // Read in the rest of configuration
  float end_time = GetXMLElementDouble({"Afterburner", "SMASH", "end_time"});
  config["General"]["End_Time"] = end_time;
  only_final_decays_ =
      GetXMLElementInt({"Afterburner", "SMASH", "only_decays"});
  JSINFO << "End time for SMASH is set to " << end_time << " fm/c";
  if (only_final_decays_) {
    JSINFO << "SMASH will only perform resonance decays, no propagation";
  }
  // output path is just dummy here, because no output from SMASH is foreseen
  JSINFO << "Seting up SMASH Experiment object";
  boost::filesystem::path output_path("./smash_output");
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
  modus->jetscape_hadrons_ = soft_particlization_sampler_->Hadron_list_;
  const int n_events = modus->jetscape_hadrons_.size();
  JSINFO << "SMASH: obtained " << n_events << " events from particlization";
  smash::Particles *smash_particles = smash_experiment_->particles();
  for (unsigned int i = 0; i < n_events; i++) {
    JSINFO << "Event " << i << " SMASH starts with "
           << modus->jetscape_hadrons_[i].size() << " particles.";
    smash_experiment_->initialize_new_event(i);
    if (!only_final_decays_) {
      smash_experiment_->run_time_evolution();
    }
    smash_experiment_->do_final_decays();
    smash_experiment_->final_output(i);
    smash_particles_to_JS_hadrons(*smash_particles,
                                  modus->jetscape_hadrons_[i]);
    JSINFO << modus->jetscape_hadrons_[i].size() << " hadrons from SMASH.";
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
    const int hadron_status = -1;
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
