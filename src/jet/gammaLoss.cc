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

#include "gammaLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"

#include <string>

#include <iostream>

#include "FluidDynamics.h"
#include <GTL/dfs.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

const double QS = 0.9;

// Register the module with the base class
RegisterJetScapeModule<gammaLoss> gammaLoss::reg("gammaLoss");

static Pythia8::Pythia PythiaFunction("IntentionallyEmpty", false);

bool gammaLoss::flag_init = 0;

gammaLoss::gammaLoss() {
  SetId("gammaLoss");
  VERBOSE(8);
  gammaLoss_on = true;
  
  iEvent = 0;
  NUM1 = 0;
}

gammaLoss::~gammaLoss() { VERBOSE(8); }

void gammaLoss::Init() {
  JSINFO << "Initialize gammaLoss ...";

  in_vac = false;
  brick_med = true;
  recoil_on = false;

  int flagInt = -100;
  double inputDouble = -99.99;

  gammaLoss_on = GetXMLElementInt({"Eloss", "gammaLoss", "gammaLoss_on"});

  JSINFO << MAGENTA << "gammaLoss input parameter";
  JSINFO << MAGENTA << "gammaLoss shower on: " << gammaLoss_on;

  //...initialize the random number generator
  srand((unsigned)time(NULL));
  NUM1 = -1 * rand();
  //    NUM1=-33;
  iEvent = 0;
}

void gammaLoss::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("ElossModule Parton List: " + GetId());
  f->WriteComment("Energy loss to be implemented accordingly ...");
}

void gammaLoss::Dump_pIn_info(int i, vector<Parton> &pIn) {
  JSWARN << "i=" << i << " gammaLoss -- status: " << pIn[i].pstat()
         << " color: " << pIn[i].color() << "  " << pIn[i].anti_color();
  JSWARN << "pid = " << pIn[i].pid() << " E = " << pIn[i].e()
         << " px = " << pIn[i].p(1) << " py = " << pIn[i].p(2)
         << "  pz = " << pIn[i].p(3) << " virtuality = " << pIn[i].t()
         << " form_time in fm = " << pIn[i].form_time()
         << " split time = " << pIn[i].form_time() + pIn[i].x_in().t();
}

void gammaLoss::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn, vector<Parton> &pOut){
  for(auto & thisparton : pIn){
    if(thisparton.pid() == 22) JSINFO << "Photon found with label " << thisparton.plabel() << " and status " << thisparton.pstat();
  }
  return;
}