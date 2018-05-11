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

#ifndef ELOSSVALIDATION_H
#define ELOSSVALIDATION_H

#include "JetEnergyLossModule.h"

using namespace Jetscape;

class ElossValidate : public JetEnergyLossModule<ElossValidate>
{  
 public:
  
  ElossValidate();
  virtual ~ElossValidate();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);


 protected:
  
};


#endif // ELOSSVALIDATION_H

