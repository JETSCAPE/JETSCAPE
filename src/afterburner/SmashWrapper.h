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

#ifndef SMASHWRAPPER_H
#define SMASHWRAPPER_H

#include "Afterburner.h"
#include "JetScapeWriter.h"

#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/listmodus.h"

using namespace Jetscape;

class SmashWrapper: public Afterburner {
 private:
    std::vector<std::vector<shared_ptr<Hadron>>> final_hadrons_;
    tinyxml2::XMLElement *smash_xml_;
    bool only_final_decays_ = false;
    shared_ptr<smash::Experiment<smash::AfterburnerModus>> smash_experiment_;

    void smash_particles_to_JS_hadrons(const smash::Particles& smash_particles,
                                  std::vector<shared_ptr<Hadron>>& JS_hadrons);

    void JS_hadrons_to_smash_particles(
        const std::vector<shared_ptr<Hadron>>& JS_hadrons,
        smash::Particles& smash_particles);
 public:
    SmashWrapper();
    void InitTask();
    void ExecuteTask();
    void WriteTask(weak_ptr<JetScapeWriter> w);
};

#endif  // SMASHWRAPPER_H
