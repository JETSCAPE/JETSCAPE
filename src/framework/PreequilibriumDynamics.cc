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

    JSINFO <<"Intialize PreequilibriumDynamics : " << GetId() << " ...";

    fd = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement(
                                                        "Preequilibrium");

    if (!fd) {
        JSWARN << "Not a valid JetScape XML Preequilibrium Dynamics section file "
           << "or no XML file loaded!";
        exit(-1);
    }

    VERBOSE(8);

    // this is grabbing the initial entropy density ?
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        JSWARN << "No initialization module, try: "
             << "auto trento = make_shared<TrentoInitial>(); "
             << "jetscape->Add(trento);";
    }

    InitializePreequilibrium(parameter_list_);

    InitTask();

    JetScapeTask::InitTasks();
}

void PreequilibriumDynamics::Exec() {
    JSINFO << "Run Preequilibrium : " << GetId() << " ...";
    VERBOSE(8) << "Current Event #" << GetCurrentEvent();

    if (ini) {
      VERBOSE(3) << "length of entropy density vector=" << ini->GetEntropyDensityDistribution().size();
    }

    EvolvePreequilibrium();

    JetScapeTask::ExecuteTasks();
}

void PreequilibriumDynamics::Clear() {
    e_.clear();
    P_.clear();
    utau_.clear();
    ux_.clear();
    uy_.clear();
    ueta_.clear();
    pi00_.clear();
    pi01_.clear();
    pi02_.clear();
    pi03_.clear();
    pi11_.clear();
    pi12_.clear();
    pi13_.clear();
    pi22_.clear();
    pi23_.clear();
    pi33_.clear();
    bulk_Pi_.clear();
}

}  // end namespace Jetscape
