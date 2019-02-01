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
// This is a causal liquefier with the JETSCAPE framework
// -----------------------------------------

#include "causal_liquefier.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

Causal_Liquefier::Causal_Liquefier() {
    SetId("Causal Liquefier");
}


void Causal_Liquefier::Init() {
    JSINFO << "Initialize causal liquefier";
}

void Causal_Liquefier::Exec() {
    JSINFO << "running causal liquefier";
}

void Causal_Liquefier::Clear() {
    JSINFO << "Finish causal liquefier";
}
