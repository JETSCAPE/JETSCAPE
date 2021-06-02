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

// Create a pythia collision at a specified point and return the two inital hard partons

#ifndef INITIALSTATERADIATIONTEST_H
#define INITIALSTATERADIATIONTEST_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class InitialStateRadiationTest : public HardProcess, public Pythia8::Pythia {

public:

  ~InitialStateRadiationTest();

  void InitTask();
  void Exec();

  virtual any GetHistory() {return any(pShowerMaster);}

private:

  const double eps = 1e-5;
  const int timeLike_stat = 23;
  const int spaceLike_stat = 25;

  // dummy variables
  unsigned int n_timeStep = 2;
  double deltaT = 0.1;
  double currentTime;

  // vector of earliest time for each parton shower
  vector<double> timeVec;

  // vector of space-like partons
  vector<shared_ptr<PartonShower>> pShowerMaster;

  node vStart, vEnd;
  vector<map<node, Parton>> nodePartonPairVec;

  const int droplet_stat = -11;
  const int miss_stat = -13;
  const int neg_stat = -17;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<InitialStateRadiationTest> reg;

  void BackwardISR();
  void ForwardISR();

};

#endif // INITIALSTATERADIATIONTEST_H
