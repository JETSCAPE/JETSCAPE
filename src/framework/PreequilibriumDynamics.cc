/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "PreequilibriumDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

PreequilibriumDynamics::PreequilibriumDynamics() {
    VERBOSE(8);
    SetId("PreequilibriumDynamics");
}

PreequilibriumDynamics::~PreequilibriumDynamics() {
    VERBOSE(8);
    disconnect_all();
}

void PreequilibriumDynamics::Init() {
    JetScapeModuleBase::Init();

    INFO <<"Intialize PreequilibriumDynamics : " << GetId() << " ...";

    fd = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement(
                                                        "Preequilibrium");

    if (!fd) {
        WARN << "Not a valid JetScape XML Preequilibrium Dynamics section file "
           << "or no XML file loaded!";
        exit(-1);
    }

    VERBOSE(8);

    // this is grabbing the initial entropy density ?
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        WARN << "No initialization module, try: "
             << "auto trento = make_shared<TrentoInitial>(); "
             << "jetscape->Add(trento);";
    }

    InitializePreequilibrium(parameter_list_);

    InitTask();

    JetScapeTask::InitTasks();
}

void PreequilibriumDynamics::Exec() {
    INFO << "Run Preequilibrium : " << GetId() << " ...";
    VERBOSE(8) << "Current Event #" << GetCurrentEvent();

    if (ini) {
      VERBOSE(3) << "length of entropy density vector=" << ini->GetEntropyDensityDistribution().size();
    }

    EvolvePreequilibrium();

    JetScapeTask::ExecuteTasks();
}

}  // end namespace Jetscape
