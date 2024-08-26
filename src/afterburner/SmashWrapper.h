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

#ifndef SMASHWRAPPER_H
#define SMASHWRAPPER_H

#include "Afterburner.h"
#include "JetScapeWriter.h"
#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/listmodus.h"

using namespace Jetscape;

/**
 * An Afterburner modus, which acts similarly to SMASH ListModus,
 * initializes from Jetscape hadrons, not from file. Needed to use
 * SMASH as a 3rd party Afterburner.
 */
class AfterburnerModus : public smash::ListModus {
 public:
  // Unlike for ListModus there is no need to get any data from the config
  AfterburnerModus(smash::Configuration config,
                   const smash::ExperimentParameters &) {
    JSINFO << "Constructing AfterburnerModus";
    config.clear();
  }
  void reset_event_numbering() { event_number_ = 0; }
  // The converter is not static, because modus holds int variables
  // for the number of warnings, which are used in try_create_particle,
  // called by this function. Maybe I (oliiny) will change this design in SMASH
  // later, but now I have to put this converter inside the AfterburnerModus.
  void JS_hadrons_to_smash_particles(
      const std::vector<shared_ptr<Hadron>> &JS_hadrons,
      smash::Particles &smash_particles);

  // This function overrides the function from ListModus.
  double initial_conditions(smash::Particles *particles,
                            const smash::ExperimentParameters &) {
    JS_hadrons_to_smash_particles(jetscape_hadrons_[event_number_], *particles);
    backpropagate_to_same_time(*particles);
    event_number_++;
    return start_time_;
  }
  std::vector<std::vector<shared_ptr<Hadron>>> jetscape_hadrons_;

 private:
  int event_number_ = 0;
};

class SmashWrapper : public Afterburner {
 private:
  bool only_final_decays_ = false;
  double end_time_ = -1.0;
  shared_ptr<smash::Experiment<AfterburnerModus>> smash_experiment_;

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<SmashWrapper> reg;

 public:
  void smash_particles_to_JS_hadrons(
      const smash::Particles &smash_particles,
      std::vector<shared_ptr<Hadron>> &JS_hadrons);
  SmashWrapper();
  void InitTask();
  void ExecuteTask();
  void WriteTask(weak_ptr<JetScapeWriter> w);
};

#endif  // SMASHWRAPPER_H
