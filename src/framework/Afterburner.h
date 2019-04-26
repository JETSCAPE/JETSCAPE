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

#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#include "JetScapeModuleBase.h"
#include "SoftParticlization.h"
#include "tinyxml2.h"
#include "sigslot.h"

namespace Jetscape {

  /// Interface to hadronic afterburner
  class Afterburner : public JetScapeModuleBase {
   public:
    Afterburner() {
      VERBOSE(8);
      SetId("Afterburner");
    }

    ~Afterburner() {
      VERBOSE(8);
      disconnect_all();
    }

    tinyxml2::XMLElement* XMLConfiguration() const { return xml_config_; }

    virtual void Init();
    virtual void Exec();
   protected:
    /// Points to options specific to the afterburner
    tinyxml2::XMLElement *xml_config_;
    /// Pointer to particlization sampler, which provides initial hadrons
    std::shared_ptr<SoftParticlization> soft_particlization_sampler_;
  };

} // end namespace Jetscape

#endif  // AFTERBURNER_H
