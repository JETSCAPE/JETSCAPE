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

#include "CausalLiquefier.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

CausalLiquefier::CausalLiquefier() {
    SetId("Causal Liquefier");
}


void CausalLiquefier::Init() {
    JSINFO << "Initialize causal liquefier";
}

void CausalLiquefier::Exec() {
    JSINFO << "running causal liquefier";
}

void CausalLiquefier::Clear() {
    JSINFO << "Finish causal liquefier";
}
