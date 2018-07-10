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

using namespace Jetscape;

class SmashWrapper: public Afterburner {
 private:
    tinyxml2::XMLElement *smash_xml_;

 public:
    SmashWrapper();
    void InitTask();
    void ExecuteTask();
    void WriteTask(weak_ptr<JetScapeWriter> w);
};

#endif  // SMASHWRAPPER_H
