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

#include "Liquefier.h"
#include "JetScapeLogger.h"


namespace Jetscape {

Liquefier::Liquefier() {
    VERBOSE(8);
    SetId("Liquefier");
}


void Liquefier::Init() {
    JSINFO <<"Intialize Liquefier : " << GetId() << " ...";
    InitTask();
}


void Liquefier::Exec() {
    JSINFO << "Run Liquefier : " << GetId() << " ...";
    ExecuteTasks();
}


void Liquefier::get_source(Jetscape::real t, Jetscape::real &j0) {
    j0 = 0.0;
    JSDEBUG<< "source to hydro = "<< j0 <<" at t= " << t;
}

}

