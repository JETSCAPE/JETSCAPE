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
//Parton Gun Test

#include "PGun.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<PGun> PGun::reg("PGun");

Pythia8::Pythia PGun::InternalHelperPythia("IntentionallyEmpty", false);

PGun::PGun() : HardProcess() {
  fixed_pT = 0;
  parID = 21;
  flag_useHybridHad = 0;
  SetId("PGun");
  VERBOSE(8);
}

PGun::~PGun() { VERBOSE(8); }

void PGun::InitTask() {
  JSDEBUG << "Initialize PGun Brick (Test) ...";
  VERBOSE(8);

  std::string s = GetXMLElementText({"Hard", "PGun", "name"});
  JSDEBUG << s << " to be initilizied ...";

  fixed_pT = GetXMLElementDouble({"Hard", "PGun", "pT"});
  JSDEBUG << s << " with fixed pT = " << fixed_pT;
  JSINFO << "Parton Gun with fixed pT = " << fixed_pT;

  parID = GetXMLElementDouble({"Hard", "PGun", "parID"});
  JSINFO << "Parton Gun with parID = " << parID;

  flag_useHybridHad = GetXMLElementInt({"Hard", "PGun", "useHybridHad"});
  JSINFO << "Use hybrid hadronization? " << flag_useHybridHad;
}

void PGun::Exec() {
  VERBOSE(2) << "Run Hard Process : " << GetId() << " ...";

  double p[4], xLoc[4];

  double pT, rapidity, phi;
  double eta_cut = 1.0;
  double tempRand;
  const double maxN = 1.0 * RAND_MAX;
  const double PI = 3.1415926;

  double ppx, ppy, ppz, pp0, mass;

  // for (int i=0;i<1;i++)
  //    {
  //      tempRand = rand()/maxN;
  //    if(tempRand < 0.25) parID = 21;
  //    else if(tempRand < 0.50) parID = 1;
  //    else if(tempRand < 0.75) parID = 2;
  //    else parID = 3;
  //    if (parID != 21) {
  //	 tempRand = rand()/maxN;
  //	 if(tempRand < 0.50) parID = -parID;
  //      }
  //     mass = 0.0;
  mass = InternalHelperPythia.particleData.m0(parID);
  //JSINFO << BOLDYELLOW << " Mass = " << mass ;
  pT = fixed_pT; //max_pT*(rand()/maxN);

  phi = 2.0 * PI * (rand() / maxN);
  rapidity = 0; //2.0*eta_cut*(rand()/maxN)-eta_cut;
  phi = 0.0;

  p[1] = pT * cos(phi);
  p[2] = pT * sin(phi);
  p[3] = sqrt(pT * pT + mass * mass) * sinh(rapidity);
  p[0] = sqrt(pT * pT + mass * mass) * cosh(rapidity);

  // Roll for a starting point
  // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  if (!ini) {
    VERBOSE(1)
        << "No initial state module, setting the starting location to 0. ";
  } else {
    double x, y;
    ini->SampleABinaryCollisionPoint(x, y);
    xLoc[1] = x;
    xLoc[2] = y;
  }

  xLoc[1] = 0.0;
  xLoc[2] = 0.0;

  if (flag_useHybridHad != 1) {
    AddParton(make_shared<Parton>(0, parID, 0, pT, rapidity, phi, p[0], xLoc));
  } else {
    auto ptn = make_shared<Parton>(0, parID, 0, pT, rapidity, phi, p[0], xLoc);
    ptn->set_color((parID > 0) ? 100 : 0);
    ptn->set_anti_color(((parID > 0) || (parID == 21)) ? 0 : 101);
    ptn->set_max_color(102);
    AddParton(ptn);
  }

  //      VERBOSEPARTON(7,*GetPartonAt(i)) <<" added "<<" at x=" << xLoc[1]<<", y=" << xLoc[2]<<", z=" << xLoc[3];
  //    }

  VERBOSE(8) << GetNHardPartons();
}
