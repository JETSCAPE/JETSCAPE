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

#ifndef NULLPREDYNAMICS_H
#define NULLPREDYNAMICS_H

#include "PreequilibriumDynamics.h"

using namespace Jetscape;

class NullPreDynamics: public PreequilibriumDynamics {
  
 public:
  NullPreDynamics();
  ~NullPreDynamics() {};
  
  void InitializePreequilibrium() {};
  void EvolvePreequilibrium();
  
 private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<NullPreDynamics> reg;
};

#endif  // NULLPREDYNAMICS_H
