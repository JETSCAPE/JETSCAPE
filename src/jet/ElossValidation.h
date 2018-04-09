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

