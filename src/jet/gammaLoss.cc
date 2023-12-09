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
#include "tinyxml2.h"
#include "RtypesCore.h"
#include "TF1.h"
#include "Math/DistSampler.h"
#include "Math/Factory.h"

#include <string>

#include <iostream>

#include "FluidDynamics.h"
#include <GTL/dfs.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

//constants
const double QS = 0.9;
const double alpha = 1./137.;
const double alphas = 0.3;

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

  integral1 = 0;
  integral2 = 0;
  x0 = 0;
  x1 = 0;
}

gammaLoss::~gammaLoss() { VERBOSE(8); }

void gammaLoss::Init() {
  JSINFO << "Initialize gammaLoss ...";

  in_vac = false;
  brick_med = true;
  recoil_on = false;
  reabsorption = true;
  hydro_Tc = 0.16;
  brick_length = 5.0;
  ratesource = 1;
  emissionOn = 0;
  iEvent = 0;

  int flagInt = -100;
  double inputDouble = -99.99;
  unsigned int xml_uintin = 0;
  int rand_seed = 0;

  gammaLoss_on = GetXMLElementInt({"Eloss", "gammaLoss", "gammaLoss_on"});
  hydro_Tc = GetXMLElementDouble({"Eloss", "Matter", "hydro_Tc"});
  brick_length = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});
  ratesource = GetXMLElementDouble({"Eloss", "gammaLoss", "source"});
  emissionOn = GetXMLElementDouble({"Eloss", "gammaLoss", "thermalEmission"});
  reabsorption = GetXMLElementDouble({"Eloss", "gammaLoss", "reabsorption"});
  x0 = GetXMLElementDouble({"Eloss", "gammaLoss", "infraredCutOff"});
  x1 = GetXMLElementDouble({"Eloss", "gammaLoss", "maxThermalEnergy"});

  JSINFO << MAGENTA << "gammaLoss input parameter";
  JSINFO << MAGENTA << "gammaLoss shower on: " << gammaLoss_on;
  JSINFO << MAGENTA << "gammaLoss rate: " << ratesource;
  JSINFO << MAGENTA << "Thermal Emission : " << emissionOn;

  //...initialize the random number generator
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  if ( RandomXmlDescription ){
    tinyxml2::XMLElement *xmle; 
    xmle = RandomXmlDescription->FirstChildElement( "seed" );
    xmle->QueryUnsignedText(&xml_uintin);
    if(xml_uintin < std::numeric_limits<unsigned int>::max()) rand_seed = xml_uintin;
  }else{
    JSWARN << "No <Random> element found in xml, seeding to 0";
    rand_seed = (unsigned)time(NULL);
  }

  //default handling
  if(xml_uintin == 0) rand_seed = (unsigned)time(NULL);

  srand(rand_seed);
  NUM1 = -1 * rand();
  rng_engine.seed(rand_seed);

  //initial integral initializations
  TF1* ffunc = new TF1("ffunc",f,x0,x1,0);
  ffunc->SetParameter(0,1);
  integral1 = ffunc->Integral(x0,x1);
  delete ffunc;

  TF1* gfunc = new TF1("gfunc",g,x0,x1,0);
  gfunc->SetParameter(0,1);
  integral2 = gfunc->Integral(x0,x1);
  delete gfunc;

  //distribution initialization
  thermalpdf = new TF1("pdf",dRdx,x0,x1,1);
  thermalpdf->SetParameter(0,0.25);
  sampler = ROOT::Math::Factory::CreateDistSampler();
  sampler->SetFunction(*thermalpdf,1);
  sampler->SetRange(x0,x1);
  sampler->SetSeed(rand_seed);
  bool success = sampler->Init();
  if(!success) JSWARN << "Sampler failed to initialize";
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
  //JSINFO << pIn.size();

  //initial declarations
  double velocity[4], xStart[4];
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  //parton loop
  for(int i=0; i<pIn.size(); i++){
    if(pIn[i].pid() != 22) continue;
    //JSINFO << "Photon found with label " << pIn[i].plabel() << " and status " << pIn[i].pstat();
    if(abs(pIn[i].pstat()) == 22) continue; //skipping absorbed photons and final state photons
    if((pIn[i].pstat() == 23 || pIn[i].pstat() == 24) && reabsorption == false){
      JSDEBUG << "Skipping thermal photon";
      continue;
    }

    //do thermal emission if triggered
    if(pIn[i].pstat() == -23){
      JSDEBUG << "Running thermal step";
      doEmission(pIn, pOut, deltaT, time);
      //JSINFO << pIn.size() << " " << pOut.size();
      continue;
    }

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
        length = 0;
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
    double deltaTprime = deltaT/tLab.T() * fmToGeVinv; //added conversion to approriate time units

    //debugging statements
    //cout << "Temp: " << now_temp << ". Abs factor: " << gammaLoss::absFactor(pLab,now_temp)*100000 << endl;
    //Dump_pIn_info(i,pIn);
    
    //removing photon if its absorbed by splitting into 2 with absorption tags
    if(gammaLoss::isAbsorbed(pLab,now_temp,deltaTprime)){
      Parton *pTemp = new Parton(0,22,-22,0.0,0.0,0.0,0.0,newpos);
      pOut.push_back(*pTemp);
      pOut.push_back(*pTemp);

      if(pIn[i].pstat() == 23 || pIn[i].pstat() == 24)JSINFO << BOLDYELLOW << "Thermal photon absorbed!";
      else JSINFO << BOLDYELLOW << "Shower photon absorbed!";
    }
  }

  return;
}

//chance for photon to be absorbed from https://arxiv.org/abs/hep-ph/9405309
double gammaLoss::absFactor1(TLorentzVector pVec, double T){
  double p = pVec.P();
  return 2.0*(5.*pi/9.)*(alpha*alphas*T*T/p)*log(0.2317*p/(alphas*T));
}

//chance for photon to be absorbed from https://arxiv.org/pdf/hep-ph/0111107.pdf
double gammaLoss::absFactor2(TLorentzVector pVec, double T){
  // Constants
  double alpha = 1.0 / 137.0;
  double prfph, x, expo, fermi, C22, Cab, Ctot, p;
  
  // Calculate alpha_s at temperature 'temp'
  double alphsT = 6.0 * pi / (27.0 * log(T / 0.022));
  double gsT = sqrt(alphsT * 4.0 * pi);
  
  // Calculate x and exponential term
  p = pVec.P();
  x = pVec.P() / T;
  expo = exp(-x);
  fermi = expo / (1.0 + expo);
  
  // Calculate prfph
  prfph = alpha * alphsT * pow(T,2) * (5.0 / 9.0) / (p);
  
  // Calculate C22 and Cab
  C22 = 0.041 / x - 0.3615 + 1.01 * exp(-1.35 * x);
  Cab = sqrt(1.5) * (0.548 / pow(x, 1.5) * log(12.28 + 1.0 / x) + 0.133 * x / sqrt(1.0 + x / 16.27));
  
  // Calculate Ctot
  Ctot = 0.5 * log(2.0 * x) + C22 + Cab;
  
  // Calculate dRd3p
  double dRd3p = prfph * fermi * (log(sqrt(3.0) / gsT) + Ctot);
  return dRd3p * 4 * pow(pi,3) / expo;
}

//seeing if photon gets absorbed, flexible for adding more rates
bool gammaLoss::isAbsorbed(TLorentzVector pVec, double T, double delTime){
  double chance;

  switch(ratesource) {
    case 1:
      chance = gammaLoss::absFactor1(pVec, T)*delTime;
      break;
    case 2:
      chance = gammaLoss::absFactor2(pVec, T)*delTime;
      break;
    default:
      chance = gammaLoss::absFactor1(pVec, T)*delTime;
      break;
  }


  if((float)rand()/RAND_MAX < chance) return true;
  else return false;
}

void gammaLoss::doEmission(vector<Parton> &pIn, vector<Parton> &pOut, double deltaT, double time){
  //initial declarations
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
  double increment = 1.0; //distance increment 
  double maxL = brick_length - (0.5*increment); //max cube distance
  double now_temp;
  int photonsmade = 0;
  double volumesampled = 0.0;

  //spacial loop
  for(double x = -1.0*maxL; x<=maxL; x+=increment){
    for(double y = -1.0*maxL; y<=maxL; y+=increment){
      for(double z = -1.0*maxL; z<=maxL; z+=increment){
        //JSINFO << "new cell at time " << time << " at position " << x << " " << y << " " << z << " ";
        GetHydroCellSignal(time, x, y, z, check_fluid_info_ptr);
        now_temp = check_fluid_info_ptr->temperature;
        //JSINFO << now_temp;
        if(now_temp < 0.1){
          //JSINFO << "Temp under threshold and skipping cell: " << now_temp;
          continue;
        }
        double posVector[4] = {time,x,y,z};

        //Lorentz math for boosting
        TLorentzVector tLab(increment,increment,increment,deltaT);
        TVector3 vMed(abs(check_fluid_info_ptr->vx), abs(check_fluid_info_ptr->vy), abs(check_fluid_info_ptr->vz));
        //TVector3 vMed(0.8, 0, 0);
        tLab.Boost(vMed); 
        //tLab.Print();
        tLab *= fmToGeVinv;
        volumesampled += tLab(0)*tLab(1)*tLab(2)*tLab(3);

        for(int i=0; i<photonsProduced(tLab, now_temp); i++){
          //JSINFO << "Making photon";
          photonsmade++;
          pIn.push_back(makeThermalPhoton(now_temp, vMed, posVector));
          JSINFO << BOLDYELLOW << "Photon made thermally!";
        }
      }
    }
  }

  //JSINFO << "Volume sampled: " << volumesampled;
  //JSINFO << "Photons made: " << photonsmade;
}

int gammaLoss::photonsProduced(TLorentzVector cell, double temp){
  double volume = cell(0)*cell(1)*cell(2);
  double deltat = cell(3);

  double gammaTot = volume*deltat*B(temp)*(log(sqrt(3)/gS(temp))*integral1 + integral2);

  //random generation samplling
  std::poisson_distribution<int> poisson_gamma(gammaTot);
  return poisson_gamma(getRandomGenerator());
}

Parton gammaLoss::makeThermalPhoton(double temp, TVector3 vMed, double position[]){
  //setting sampled quantities
  double momentum = temp*sampler->Sample1D();
  TLorentzVector partmom(momentum,0,0,momentum);
  partmom.SetTheta(genTheta());
  partmom.SetPhi(genPhi());

  //boosting to lab frame then output
  partmom.Boost(vMed); 
  FourVector jsP(partmom.Px(),partmom.Py(),partmom.Pz(),partmom.E());
  Parton *pTemp = new Parton(0,22,23,jsP,position);
  return *pTemp;
}

double gammaLoss::alphaS(double temp){
  return 6.0 * pi / (27.0 * log(temp / 0.022));
}

double gammaLoss::gS(double temp){
  return sqrt(4.0 * pi * alphaS(temp));
}

//factor of T out front off emission function
double gammaLoss::B(double temp){
  return 20.*alpha*alphaS(temp)*pow(temp,4)/(9.*pi);
}

//total momentum distribtion functions
// x[0] = k/T; par[0] = T
Double_t gammaLoss::dRdx(Double_t *x, Double_t *par){
  return B(par[0])*(log(sqrt(3)/gS(par[0]))*f(x,par) + g(x,par));
}

Double_t gammaLoss::f(Double_t *x, Double_t *par){
  return x[0]/(1.0+exp(x[0]));
}

Double_t gammaLoss::g(Double_t *x, Double_t *par){
  double C22 = (0.041/x[0]) - 0.3615 + (1.01*exp(-1.35*x[0]));
  double Cab = sqrt(1.5) * ((0.548/pow(x[0], 1.5) * log(12.28 + 1.0/x[0])) + (0.133*x[0]/sqrt(1.0 + x[0]/16.27))) ;
  return (C22 + Cab + 0.5*log(2.0*x[0]))*f(x,par);
}