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

#ifndef KTISOWRAPPER_H
#define KTISOWRAPPER_H

#include "PreequilibriumDynamics.h"
#include "KTiso.cpp" 

using namespace Jetscape;

class KTIsoWrapper: public PreequilibriumDynamics {
 private:
  //int mode; //!< records running mode
  ITA *ktiso_ptr;
  
 public:
  KTIsoWrapper();
  ~KTIsoWrapper();
  
  void InitializePreequilibrium(PreEquilibriumParameterFile parameter_list);  
  void EvolvePreequilibrium();
};

#endif  // KTISOWRAPPER_H
