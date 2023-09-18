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
}

void gammaLoss::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn, vector<Parton> &pOut){
  //JSINFO << "gamma";

  //initial declarations
  double velocity[4], xStart[4], velocity_jet[4];
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  //parton loop
  for(int i=0; i<pIn.size(); i++){
    if(pIn[i].pid() == 22) JSINFO << "Photon found with label " << pIn[i].plabel() << " and status " << pIn[i].pstat();
    if(pIn[i].pid() != 22) continue;

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

    // Define a vector in the direction of the jet originating parton.
    // there is some amount of redundancy here, pIn[i].jet_v() is basically the jet velocity
    // we are defining a local copy in velocity_jet
    velocity_jet[0] = 1.0;
    velocity_jet[1] = pIn[i].jet_v().x();
    velocity_jet[2] = pIn[i].jet_v().y();
    velocity_jet[3] = pIn[i].jet_v().z();

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

    initRdotV = (initRx * pIn[i].jet_v().x() + initRy * pIn[i].jet_v().y() +
                 initRz * pIn[i].jet_v().z()) /
                mod_jet_v;
    initVdotV = (initVx * pIn[i].jet_v().x() + initVy * pIn[i].jet_v().y() +
                 initVz * pIn[i].jet_v().z()) /
                mod_jet_v;
    // Note: jet_v()/mod_jet_v is a unit 3 vector in the direction of the jet originating parton.

    double now_R0 = time;
    double now_Rx = initRx + (time - initR0) * initVx;
    double now_Ry = initRy + (time - initR0) * initVy;
    double now_Rz = initRz + (time - initR0) * initVz;
    double now_temp;

    double SpatialRapidity = 0.5 * std::log((now_R0 + now_Rz) / (now_R0 - now_Rz));

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

    cout << "Temp: " << now_temp << endl;
    pIn[i].set_stat(22);
    /*pOut.push_back(pIn[i]);*/
    //pIn.erase(pIn.begin()+i);
  }

  return;
}

double gammaLoss::fillQhatTab(double y) {

  double xLoc, yLoc, zLoc, tLoc;
  double vxLoc, vyLoc, vzLoc, gammaLoc, betaLoc;
  double edLoc, sdLoc;
  double tempLoc;
  double flowFactor, qhatLoc;
  int hydro_ctl;
  double lastLength = initR0;

  double tStep = 0.1;

  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  for (int i = 0; i < dimQhatTab; i++) {
    tLoc = tStep * i;

    //if(tLoc<initR0-tStep) { // potential problem of making t^2<z^2

    double boostedTStart = tStart * std::cosh(y);
    if (tLoc < initR0 || tLoc < boostedTStart) {
      qhatTab1D[i] = 0.0;
      continue;
    }

    xLoc = initRx + (tLoc - initR0) * initVx;
    yLoc = initRy + (tLoc - initR0) * initVy;
    zLoc = initRz + (tLoc - initR0) * initVz;

    //        if(bulkFlag == 1) { // read OSU hydro
    //            readhydroinfoshanshan_(&tLoc,&xLoc,&yLoc,&zLoc,&edLoc,&sdLoc,&tempLoc,&vxLoc,&vyLoc,&vzLoc,&hydro_ctl);
    //        } else if(bulkFlag == 2) { // read CCNU hydro
    //            hydroinfoccnu_(&tLoc, &xLoc, &yLoc, &zLoc, &tempLoc, &vxLoc, &vyLoc, &vzLoc, &hydro_ctl);
    //        } else if(bulkFlag == 0) { // static medium
    //            vxLoc = 0.0;
    //            vyLoc = 0.0;
    //            vzLoc = 0.0;
    //            hydro_ctl = 0;
    //            tempLoc = T;
    //        }

    if (std::isinf(tLoc) || std::isnan(tLoc) || std::isinf(zLoc) ||
        std::isnan(zLoc) || std::abs(zLoc) > tLoc) {
      JSWARN << "Third instance";
      JSWARN << "Loc for vector is:" << tLoc << ", " << xLoc << ", " << yLoc
             << ", " << zLoc;
      JSWARN << "initR0, initRx, initRy, initRz="
             << ", " << initR0 << ", " << initRx << ", " << initRy << ", "
             << initRz;
      JSWARN << "initVx, initVy, initVz =" << initVx << ", " << initVy << ", "
             << initVz;
      JSWARN << "initVMod=" << std::setprecision(20)
             << std::sqrt(initVx * initVx + initVy * initVy + initVz * initVz);
      JSWARN << "Can't dump pIn_info as we are in fillQhatTab. But it should "
                "be dumped right before this."; //Dump_pIn_info(i, pIn);
                                                //exit(0);
    }

    GetHydroCellSignal(tLoc, xLoc, yLoc, zLoc, check_fluid_info_ptr);
    VERBOSE(8) << MAGENTA << "Temperature from medium = "
               << check_fluid_info_ptr->temperature;

    tempLoc = check_fluid_info_ptr->temperature;
    sdLoc = check_fluid_info_ptr->entropy_density;
    vxLoc = check_fluid_info_ptr->vx;
    vyLoc = check_fluid_info_ptr->vy;
    vzLoc = check_fluid_info_ptr->vz;

    hydro_ctl = 0;

    if (hydro_ctl == 0 && tempLoc >= hydro_Tc) {
      lastLength = tLoc;
      betaLoc = sqrt(vxLoc * vxLoc + vyLoc * vyLoc + vzLoc * vzLoc);
      gammaLoc = 1.0 / sqrt(1.0 - betaLoc * betaLoc);
      flowFactor =
          gammaLoc * (1.0 - (initVx * vxLoc + initVy * vyLoc + initVz * vzLoc));

      //if(run_alphas==1){ alphas= 4*pi/(9.0*log(2*initEner*tempLoc/0.04));}

      /* if (qhat0 < 0.0) {
        // calculate qhat with alphas
        double muD2 = 6.0 * pi * alphas * tempLoc * tempLoc;
        // if(initEner > pi*tempLoc) qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*initEner*tempLoc/muD2);
        // else qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*pi*tempLoc*tempLoc/muD2);
        // fitted formula from https://arxiv.org/pdf/1503.03313.pdf
        if (initEner > pi * tempLoc)
          qhatLoc = Ca * 50.4864 / pi * pow(alphas, 2) * pow(tempLoc, 3) *
                    log(5.7 * initEner * tempLoc / 4.0 / muD2);
        else
          qhatLoc = Ca * 50.4864 / pi * pow(alphas, 2) * pow(tempLoc, 3) *
                    log(5.7 * pi * tempLoc * tempLoc / 4.0 / muD2);
        qhatLoc = qhatLoc * flowFactor;
        if (qhatLoc < 0.0)
          qhatLoc = 0.0;
      } else { // use input qhat
        if (brick_med) {
          qhatLoc = qhat0 * 0.1973 * flowFactor;
        } else {
          qhatLoc = qhat0 / 96.0 * sdLoc * 0.1973 *
                    flowFactor; // qhat0 at s = 96fm^-3
        }
      }*/

      // GeneralQhatFunction(int QhatParametrizationType, double Temperature, double EntropyDensity, double FixAlphas,  double Qhat0, double E, double muSquare);
      double muSquare=-1;//For virtuality dependent cases, we explicitly modify q-hat inside Sudakov, due to which we set here scale=-1; Alternatively one could extend the dimension of the q-hat table

      qhatLoc= GeneralQhatFunction(QhatParametrizationType, tempLoc, sdLoc, alphas, qhat0, initEner, muSquare);      
      qhatLoc = qhatLoc * flowFactor;

      //JSINFO << "check qhat --  ener, T, qhat: " << initEner << " , " << tempLoc << " , " << qhatLoc;
    } else { // outside the QGP medium
      qhatLoc = 0.0;
    }

    qhatTab1D[i] =
        qhatLoc / sqrt(2.0); // store qhat value in light cone coordinate
  }

  for (int i = 0; i < dimQhatTab; i++) { // dim of loc

    double totValue = 0.0;

    for (int j = 0; i + j < dimQhatTab; j++) { // dim of tau_f

      totValue = totValue + qhatTab1D[i + j];
      qhatTab2D[i][j] = totValue / (j + 1);
    }
  }

  //return(lastLength*sqrt(2.0)*5.0); // light cone + GeV unit
  return ((2.0 * lastLength + initRdotV - initR0) / sqrt(2.0) *
          5.0); // light cone + GeV unit
}

//////////////////////////////////General Function of q-hat//////////////////////////////////
// E is the energy and muSquare is the virtuality of the parton
double gammaLoss::GeneralQhatFunction(int QhatParametrization, double Temperature, double EntropyDensity, double FixAlphas, double Qhat0, double E, double muSquare){
  int ActiveFlavor=3; qhat=0.0;
  double DebyeMassSquare = FixAlphas*4*pi*pow(Temperature,2.0)*(6.0 + ActiveFlavor)/6.0; 
  double ScaleNet=2*E*Temperature; 
  if(ScaleNet < 1.0){ ScaleNet=1.0; }
  switch(QhatParametrization)
    {
      //HTL formula with all alpha_s as constant and controlled by XML
    case 0:
      qhat = (Ca*50.4864/pi)*pow(FixAlphas,2)*pow(Temperature,3)*log(ScaleNet/DebyeMassSquare); 
      break;
      
      //alpha_s at scale muS=2ET and second alpha_s at muS=DebyeMassSquare is fit parameter
    case 1:	
      qhat = (Ca*50.4864/pi)*RunningAlphaS(ScaleNet)*FixAlphas*pow(Temperature,3)*log(ScaleNet/DebyeMassSquare);
      break;
      
      //Constant q-hat 
    case 2:
    qhat = Qhat0*0.1973;
    break;

    //Scale with T^3
    case 3:
    qhat = Qhat0*pow(Temperature/0.3,3)*0.1973; // w.r.t T=0.3 GeV
    break;
  
    //Scale with entropy density
    case 4:
    qhat = Qhat0*(EntropyDensity/96.0)*0.1973; // w.r.t S0=96 fm^-3
    break;
    
    //HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat  
    //Function is 1/(1+A*pow(log(Q^2),2)+B*pow(log(Q^2),4)) 
    case 5:
      qhat = (Ca*50.4864/pi)*RunningAlphaS(ScaleNet)*FixAlphas*pow(Temperature,3)*log(ScaleNet/DebyeMassSquare);
      qhat = qhat*VirtualityQhatFunction(5, E, muSquare);
    break;
    
    //HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat
    //Function is int^{1}_{xB} e^{-ax} / (1+A*pow(log(Q^2),1)+B*pow(log(Q^2),2)) 
    case 6:
      qhat = (Ca*50.4864/pi)*RunningAlphaS(ScaleNet)*FixAlphas*pow(Temperature,3)*log(ScaleNet/DebyeMassSquare);
      qhat = qhat*VirtualityQhatFunction(6, E, muSquare);
      break;

      //HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat
      //Function is int^{1}_{xB} x^{a}(1-x)^{b} / (1+A*pow(log(Q^2),1)+B*pow(log(Q^2),2)) 
    case 7:
      qhat = (Ca*50.4864/pi)*RunningAlphaS(ScaleNet)*FixAlphas*pow(Temperature,3)*log(ScaleNet/DebyeMassSquare);
      qhat = qhat*VirtualityQhatFunction(7, E, muSquare);
      break;

    default:      
      JSINFO<<"q-hat Parametrization "<<QhatParametrization<<" is not used, qhat will be set to zero";    
    }  
  return qhat;
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

//////////////////////////////////////////////////////////////////////////////////////
///////////// Virtuality dependent prefactor for q-hat function /////////////////////////
double gammaLoss::VirtualityQhatFunction(int QhatParametrization,  double enerLoc, double muSquare){
  double ans=0;double xB=0, xB0=0, IntegralNorm=0;

  if( muSquare <= Q00*Q00) {ans=1;}
  else
    {  
      switch(QhatParametrization)
	{
	case 0:
          ans =  1.0;
          break;
	  
	case 1:
          ans =  1.0;
          break;

	case 2:
          ans =  1.0;
          break;  
	  
	case 3:
          ans =  1.0;
          break;

	case 4:
          ans =  1.0;
          break;  

	case 5:
	  ans =  1.0 + qhatA*log(Q00*Q00)*log(Q00*Q00) + qhatB*pow(log(Q00*Q00),4);
	  ans = ans/( 1.0 + qhatA*log(muSquare)*log(muSquare) + qhatB*pow(log(muSquare),4)  );
	  break;
	  
	case 6: 
	  xB  = muSquare/(2.0*enerLoc);
	  xB0 = Q00*Q00/(2.0*enerLoc); if(xB<=xB0){ans=1.0; break; }
	  if ( qhatC > 0.0 && xB < 0.99)
	    {
	      ans = ( exp(qhatC*(1.0-xB)) - 1.0 )/(1.0 + qhatA*log(muSquare/0.04) + qhatB*log(muSquare/0.04)*log(muSquare/0.04) );
	      ans = ans*(1.0 + qhatA*log(Q00*Q00/0.04) + qhatB*log(Q00*Q00/0.04)*log(Q00*Q00/0.04) )/( exp(qhatC*(1.0-xB0)) - 1.0 );
	      //JSINFO<<"K xB="<<xB<<", and (E,muSquare)=("<<enerLoc<<","<<muSquare<<"), and Virtuality dep Qhat="<<ans;
	    }
	  else if( qhatC == 0.0 && xB < 0.99)
	    {
	      ans = (1.0-xB)/(1.0 + qhatA*log(muSquare/0.04) + qhatB*log(muSquare/0.04)*log(muSquare/0.04) );
	      ans = ans*(1.0 + qhatA*log(Q00*Q00/0.04) + qhatB*log(Q00*Q00/0.04)*log(Q00*Q00/0.04) )/(1-xB0);
	    }
	  else {ans=0.0;}	  
	  //JSINFO<<"L xB="<<xB<<", and (E,muSquare)=("<<enerLoc<<","<<muSquare<<"), and Virtuality dep Qhat="<<ans;
          break;

	case 7:
	  xB  = muSquare/(2.0*enerLoc);
	  xB0 = Q00*Q00/(2.0*enerLoc); if(xB<=xB0){ans=1.0; break; }
	  ans = IntegralPDF(xB, qhatC, qhatD)/(1.0 + qhatA*log(muSquare/0.04) + qhatB*log(muSquare/0.04)*log(muSquare/0.04) );
	  IntegralNorm = IntegralPDF(xB0, qhatC, qhatD);
	  if( IntegralNorm > 0.0 )
	    {
	      ans=ans*(1.0 + qhatA*log(muSquare/0.04) + qhatB*log(muSquare/0.04)*log(muSquare/0.04) )/IntegralNorm;
	    }
	  else {ans=0;}
	  //JSINFO<<"L xB="<<xB<<", and (E,muSquare)=("<<enerLoc<<","<<muSquare<<"), and Virtuality dep Qhat="<<ans; 
	  break;
	  
	default:
	  JSINFO<<"q-hat Parametrization "<<QhatParametrization<<" is not used, VirtualityQhatFunction is set to zero or one";
	}
    }

  //JSINFO<<"Qhat Type="<<QhatParametrization<<", and (E,muSquare)=("<<enerLoc<<","<<muSquare<<"), and Virtuality dep Qhat="<<ans;
  return ans;  
}

//x integration  of QGP-PDF from xB to 1
double gammaLoss::IntegralPDF(double xB, double a, double b){
  double Xmin=0.01, Xmax=0.99, X=0, dX=0, ans=0; int N=100, ix=0;
  dX = (1.0-0.0)/N;
  if(xB > Xmax) {ans=0;}
  else
    {
      ix = xB/dX;
      for(int i=ix; i<N; i++)
	{
	  X= (i+0.5)*dX;
	  if (X<Xmin){X=Xmin;}
	  ans = ans + pow(X,a)*pow(1-X,b)*dX;
	}
    }
  //JSINFO<<"(xB,a,b)=("<<xB<<","<<a<<","<<b<<"), \t Area="<<ans;
  return ans;
}