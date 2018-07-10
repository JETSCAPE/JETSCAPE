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

#include "smash/configuration.h"
#include "smash/decaymodes.h"
#include "smash/experiment.h"
#include "smash/inputfunctions.h"
#include "smash/particles.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace Jetscape;

SmashWrapper::SmashWrapper() {
  SetId("SMASH");
}

void SmashWrapper::InitTask() {
  JSINFO << "SMASH: picking SMASH-specific configuration from xml file";
  smash_xml_ = xml_config_->FirstChildElement("SMASH");
  if (!smash_xml_) {
    JSWARN << "No XML section for SMASH! Please check the input file.";
    exit(-1);
  }
  // Read SMASH config file
  string smash_config = smash_xml_->
      FirstChildElement("SMASH_config_file")->GetText();
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
  string hadron_list = smash_xml_->
      FirstChildElement("SMASH_particles_file")->GetText();
  if (boost::filesystem::exists(hadron_list)) {
    config["particles"] =
        smash::read_all(boost::filesystem::ifstream{hadron_list});
    smash::ParticleType::create_type_list(config.take({"particles"}));
    JSINFO << "Obtaining SMASH hadron list from " << hadron_list;
  } else {
    JSINFO << "Using default SMASH hadron list.";
  }
  // SMASH decaymodes
  string decays_list = smash_xml_->
      FirstChildElement("SMASH_decaymodes_file")->GetText();
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
  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config["General"]["Randomseed"] = random_seed;
  VERBOSE(2) << "Random seed used for SMASH: " << random_seed;
  // Set SMASH logging
  smash::set_default_loglevel(
      config.take({"Logging", "default"}, einhard::INFO));
  smash::create_all_loggers(config["Logging"]);
}

void SmashWrapper::ExecuteTask() {
  // todo(oliiny): Is it really avoiding copy?
  const auto &initial_hadrons = soft_particlization_sampler_->Hadron_list_;
  const int n_events = initial_hadrons.size();
  JSINFO << "SMASH: obtained " << n_events << " events from particlization";
  for (unsigned int i = 0; i < n_events; i++) {
    JSINFO << initial_hadrons[i].size() << " particles in event " << i;
  }
}

void SmashWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSINFO << "SMASH hadronic afterburner printout";
}
