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

//photons given code 22 sent to final state
//photons given code -22 removed entirely

#include "gammaLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"
#include "Math/Boost.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <string>

#include <iostream>

#include "FluidDynamics.h"
#include <GTL/dfs.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

const double QS = 0.9;
const double alpha = 1./137.;
const double alphaS = 0.3;

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
  hydro_Tc = 0.;
  brick_length = 0.;
  initR0 = 0.;
  initRx = 0.;
  initRy = 0.;
  initRz = 0.;
  initVx = 0.;
  initVy = 0.;
  initVz = 0.;
  initRdotV = 0.;
  initVdotV = 0.;
  initEner = 0.;
  qhat0 = 0.;
  alphas = 0.;
  tscale=1;
  QhatParametrizationType=-1;
  qhatA=0.;
  qhatB=0.;
  qhatC=0.;
  qhatD=0.;
}

gammaLoss::~gammaLoss() { VERBOSE(8); }

void gammaLoss::Init() {
  JSINFO << "Initialize gammaLoss ...";

  in_vac = false;
  brick_med = true;
  recoil_on = false;
  hydro_Tc = 0.16;
  brick_length = 4.0;
  qhat = 0.0;
  Q00 = 1.0;    // virtuality separation scale
  qhat0 = 2.0;  // GeV^2/fm for gluon at s = 96 fm^-3
  alphas = 0.3; // only useful when qhat0 is a negative number
  tscale=1;
  QhatParametrizationType=-1;
  qhatA=1;
  qhatB=1;
  qhatC=1;
  qhatD=1;

  int flagInt = -100;
  double inputDouble = -99.99;

  gammaLoss_on = GetXMLElementInt({"Eloss", "gammaLoss", "gammaLoss_on"});
  hydro_Tc = GetXMLElementDouble({"Eloss", "Matter", "hydro_Tc"});
  brick_length = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});

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
  JSWARN << "Position: " << pIn[i].x_in().x() << " " << pIn[i].x_in().y()  << " "<< pIn[i].x_in().z() ;
}

void gammaLoss::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn, vector<Parton> &pOut){
  //JSINFO << "gamma";

  //initial declarations
  double velocity[4], xStart[4];
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  //parton loop
  for(int i=0; i<pIn.size(); i++){
    if(pIn[i].pid() != 22) continue;
    //JSINFO << "Photon found with label " << pIn[i].plabel() << " and status " << pIn[i].pstat();
    if(abs(pIn[i].pstat()) == 22) continue; //skipping absorbed photons and final state photons

    //velocity and spatial settings
    velocity[0] = 1.0;
    // Define 3 velocity of the parton
    for (int j = 1; j <= 3; j++) {
      velocity[j] = pIn[i].p(j) / pIn[i].e();
    }
    // velocityMod will be 1 for most partons except maybe heavy quarks, we say "maybe", as for very
    // energetic heavy quarks, the velocity may be very close to 1.
    double velocityMod =
        std::sqrt(std::pow(velocity[1], 2) + std::pow(velocity[2], 2) +
                  std::pow(velocity[3], 2));

    if (velocityMod > 1.0 + rounding_error) {
      JSINFO << BOLDRED
             << " tachyonic propagation detected for parton passed from hard "
                "scattering, velocity mod = "
             << velocityMod;
      JSWARN << "velocityMod=" << std::setprecision(20) << velocityMod;
      Dump_pIn_info(i, pIn);
      //assert(velocityMod < 1.0 + rounding_error);
    }

    if (pIn[i].form_time() < 0.0) pIn[i].set_jet_v(velocity); // jet velocity is set only once
    // Notice the assumption that partons passed from hard scattering are on shell.
    // If they had virtuality then this would not be correct.

    // Modulus of the vector pIn[i].jet_v
    double mod_jet_v =
        std::sqrt(pow(pIn[i].jet_v().x(), 2) + pow(pIn[i].jet_v().y(), 2) +
                  pow(pIn[i].jet_v().z(), 2));

    for (int j = 0; j <= 3; j++) {
      xStart[j] = pIn[i].x_in().comp(j);
    }

    // SC: read in hydro
    initR0 = xStart[0];
    initRx = xStart[1];
    initRy = xStart[2];
    initRz = xStart[3];
    initVx = velocity[1] / velocityMod;
    initVy = velocity[2] / velocityMod;
    initVz = velocity[3] / velocityMod;

    if (std::abs(pIn[i].pid()) == 4 || std::abs(pIn[i].pid()) == 5) {
      double OnShellEnergy = std::sqrt(
          pIn[i].px() * pIn[i].px() + pIn[i].py() * pIn[i].py() +
          pIn[i].pz() * pIn[i].pz() + pIn[i].restmass() * pIn[i].restmass());

      initVx = pIn[i].px() / OnShellEnergy;
      initVy = pIn[i].py() / OnShellEnergy;
      initVz = pIn[i].pz() / OnShellEnergy;

      velocityMod =
          std::sqrt(initVx * initVx + initVy * initVy + initVz * initVz);
    }

    //current position
    double now_R0 = time;
    double now_Rx = initRx + (time - initR0) * initVx;
    double now_Ry = initRy + (time - initR0) * initVy;
    double now_Rz = initRz + (time - initR0) * initVz;
    double now_temp;

    double SpatialRapidity = 0.5 * std::log((now_R0 + now_Rz) / (now_R0 - now_Rz));
    
    //updating pos info
    double newpos[4] = {now_R0,now_Rx,now_Ry,now_Rz};
    pIn[i].set_x(newpos);

    //hydro settings
    double length;
    double initEner = pIn[i].e(); // initial Energy of parton
    if (!in_vac) {
      if (GetJetSignalConnected())
        length = 0;//fillQhatTab(SpatialRapidity);
      else {
        JSWARN << "Couldn't find a hydro module attached!";
        throw std::runtime_error(
            "Please attach a hydro module or set in_vac to 1 in the XML file");
      }
    }
    if(brick_med) length = brick_length*fmToGeVinv; /// length in GeV-1 will have to changed for hydro

    //getting temp from hydro
    double boostedTStart = tStart * cosh(SpatialRapidity);
    if (!in_vac && now_R0 >= boostedTStart) {
      if (now_R0 * now_R0 < now_Rz * now_Rz)
        cout << "Warning 1: " << now_R0 << "  " << now_Rz << endl;
      GetHydroCellSignal(now_R0, now_Rx, now_Ry, now_Rz, check_fluid_info_ptr);
      //VERBOSE(8)<<MAGENTA<<"Temperature from medium = "<<check_fluid_info_ptr->temperature;
      now_temp = check_fluid_info_ptr->temperature;
      //JSINFO << BOLDYELLOW << "MATTER time = " << now_R0 << " x = " << now_Rx << " y = " << now_Ry << " z = " << now_Rz << " temp = " << now_temp;
      //JSINFO << BOLDYELLOW << "MATTER initVx, initVy, initVz =" << initVx << ", " << initVy << ", " << initVz;
      //JSINFO << BOLDYELLOW << "MATTER velocityMod=" << velocityMod;
    } else {
      now_temp = 0.0;
    }

    //sending photons that have left medium to the final state
    if(now_temp == 0.0){
      pIn[i].set_stat(22);
      continue;
    }

    //Lorentz math for boosting
    TLorentzVector pLab(pIn[i].px(),pIn[i].py(),pIn[i].pz(),pIn[i].e());
    TLorentzVector tLab(0.,0.,0.,1.);
    TVector3 vMed(check_fluid_info_ptr->vx, check_fluid_info_ptr->vy, check_fluid_info_ptr->vz);
    pLab.Boost(-vMed);
    tLab.Boost(-vMed); 
    double deltaTprime = deltaT/tLab.T();

    //debugging statements
    //cout << "Temp: " << now_temp << ". Abs factor: " << gammaLoss::absFactor(pLab,now_temp)*100000 << endl;
    //Dump_pIn_info(i,pIn);
    
    //removing photon if its absorbed
    if(gammaLoss::isAbsorbed(pLab,now_temp,deltaTprime)){
      Parton *pTemp = new Parton(0,22,-22,0.0,0.0,0.0,0.0,newpos);
      pOut.push_back(*pTemp);
      pOut.push_back(*pTemp);
      JSINFO << BOLDYELLOW << "Photon absorbed!";
    }
  }

  return;
}

//chance for photon to be absorbed from https://arxiv.org/abs/hep-ph/9405309
double gammaLoss::absFactor(TLorentzVector pVec, double T){
  double p = pVec.P();
  return 2*(5.*pi/9.)*(alpha*alphaS*T*T/p)*log(0.2317*p/(alphaS*T));
}

//seeing if photon gets absorbed
bool gammaLoss::isAbsorbed(TLorentzVector pVec, double T, double delTime){
  double chance = gammaLoss::absFactor(pVec, T)*delTime;

  if((float)rand()/RAND_MAX < chance) return true;
  else return false;
}

/////////////////// Running alphas for HTL-qhat: Do not use for others///////////////////
double gammaLoss::RunningAlphaS(double muSquare){
  int ActiveFlavor=3;
  double Square_Lambda_QCD_HTL = exp( -12.0*pi/( (33 - 2*ActiveFlavor)*alphas) );
  double ans = 12.0*pi/( (33.0- 2.0*ActiveFlavor)*log(muSquare/Square_Lambda_QCD_HTL) );
  if(muSquare < 1.0) {ans=alphas; }
  
  VERBOSE(8)<<"Fixed-alphaS="<<alphas<<", Lambda_QCD_HTL="<<sqrt(Square_Lambda_QCD_HTL)<<", mu2="<<muSquare<<", Running alpha_s"<<ans;
  return ans;
}