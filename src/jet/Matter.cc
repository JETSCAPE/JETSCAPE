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

#include "Matter.h"
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"
#include "SplittingFunction.h"
#include "SudakovFunction.h"

#include <string>

#include <iostream>

#include "FluidDynamics.h"
#include <GTL/dfs.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;
 
const double QS = 0.9;//1.0;

// Register the module with the base class
RegisterJetScapeModule<Matter> Matter::reg("Matter");

static Pythia8::Pythia PythiaFunction("IntentionallyEmpty", false);

bool Matter::flag_init = 0;

double Matter::RHQ[60][20] = {{0.0}};    //total scattering rate for heavy quark
double Matter::RHQ11[60][20] = {{0.0}};  //Qq->Qq
double Matter::RHQ12[60][20] = {{0.0}};  //Qg->Qg
double Matter::qhatHQ[60][20] = {{0.0}}; //qhat of heavy quark

double Matter::distFncB[N_T][N_p1][N_e2] = {{{0.0}}};
double Matter::distFncF[N_T][N_p1][N_e2] = {{{0.0}}};
double Matter::distMaxB[N_T][N_p1][N_e2] = {{{0.0}}};
double Matter::distMaxF[N_T][N_p1][N_e2] = {{{0.0}}};
double Matter::distFncBM[N_T][N_p1] = {{0.0}};
double Matter::distFncFM[N_T][N_p1] = {{0.0}};

Matter::Matter() {
  SetId("Matter");
  VERBOSE(8);
  qhat = 0.0;
  ehat = 0.0;
  e2hat = 0.0;
  length = 0.0;
  MaxColor = 0;
  matter_on = true;
  in_vac = false;
  brick_med = true;
  recoil_on = false;
  broadening_on = false;
  hydro_Tc = 0.;
  qhat0 = 0.;
  alphas = 0.;
  tscale=1;
  QhatParametrizationType=-1;
  qhatA=0.;
  qhatB=0.;
  qhatC=0.;
  qhatD=0.;
  brick_length = 0.;
  vir_factor = 0.;
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
  Q00 = 0.;
  Q0 = 0.;
  T0 = 0.;
  iEvent = 0;
  NUM1 = 0;
}

Matter::~Matter() { VERBOSE(8); }

void Matter::Init() {
  JSINFO << "Initialize Matter ...";

  elasticCollision.setter(0.2,0.2,0.1,55.0,0.2,10.0);

  in_vac = false;
  brick_med = true;
  recoil_on = false;

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
  hydro_Tc = 0.16;
  brick_length = 4.0;
  vir_factor = 1.0;

  double m_qhat = GetXMLElementDouble({"Eloss", "Matter", "qhat0"});
  SetQhat(m_qhat);
  //qhat = GetQhat()/fmToGeVinv ;
  //qhat0 = GetQhat()/fmToGeVinv ;
  qhat0 = GetQhat();

  int flagInt = -100;
  double inputDouble = -99.99;

  matter_on = GetXMLElementInt({"Eloss", "Matter", "matter_on"});
  in_vac = GetXMLElementInt({"Eloss", "Matter", "in_vac"});
  recoil_on = GetXMLElementInt({"Eloss", "Matter", "recoil_on"});
  broadening_on = GetXMLElementInt({"Eloss", "Matter", "broadening_on"});
  brick_med = GetXMLElementInt({"Eloss", "Matter", "brick_med"});
  Q00 = GetXMLElementDouble({"Eloss", "Matter", "Q0"});
  T0 = GetXMLElementDouble({"Eloss", "Matter", "T0"});
  alphas = GetXMLElementDouble({"Eloss", "Matter", "alphas"});
  qhatA = GetXMLElementDouble({"Eloss", "Matter", "qhatA"});
  qhatB = GetXMLElementDouble({"Eloss", "Matter", "qhatB"});
  qhatC = GetXMLElementDouble({"Eloss", "Matter", "qhatC"});
  qhatD = GetXMLElementDouble({"Eloss", "Matter", "qhatD"});
  tStart = GetXMLElementDouble({"Eloss", "tStart"});
  //run_alphas = GetXMLElementInt({"Eloss", "Matter", "run_alphas"});
  QhatParametrizationType=GetXMLElementInt({"Eloss", "Matter", "QhatParametrizationType"});
  hydro_Tc = GetXMLElementDouble({"Eloss", "Matter", "hydro_Tc"});
  brick_length = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});
  vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

  if (vir_factor < 0.0) {
    cout << "Reminder: negative vir_factor is set, initial energy will be used "
            "as initial t_max"
         << endl;
  }

  MaxColor = 101; // MK:recomb

  JSINFO << MAGENTA << "Matter input parameter";
  JSINFO << MAGENTA << "Matter shower on: " << matter_on;
  JSINFO << MAGENTA << "in_vac: " << in_vac << "  brick_med: " << brick_med
         << "  recoil_on: " << recoil_on<<", tStart ="<<tStart;
  JSINFO << MAGENTA << "Q0: " << Q00 << " vir_factor: " << vir_factor
         << "  qhat0: " << qhat0 << " alphas: " << alphas << ", QhatParametrizationType="<<QhatParametrizationType
         << "  qhatA: " << qhatA << " qhatB:  " <<qhatB  << "  qhatC: " << qhatC << " qhatD:  " <<qhatD
         << " hydro_Tc: " << hydro_Tc << " brick_length: " << brick_length;
  if(QhatParametrizationType ==1){ JSINFO <<"Running alphas will be used as 4*pi/(9.0*log(2*enerLoc*tempLoc/LambdaQCDS))"; }

  if (recoil_on && !flag_init) {
    JSINFO << MAGENTA
           << "Reminder: download LBT tables first and cmake .. if recoil is "
              "switched on in Matter.";
    read_tables(); // initialize various tables
    flag_init = true;
  }

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};

  //...initialize the random number generator
  //srand((unsigned)time(NULL));
  NUM1 =-1*static_cast<int>(ZeroOneDistribution(*GetMt19937Generator())*RAND_MAX);	
  //-1 * rand();
  //    NUM1=-33;
  iEvent = 0;
  
}

void Matter::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("ElossModule Parton List: " + GetId());
  f->WriteComment("Energy loss to be implemented accordingly ...");
}

void Matter::Dump_pIn_info(int i, vector<Parton> &pIn) {
  JSWARN << "i=" << i << " Matter -- status: " << pIn[i].pstat()
         << " color: " << pIn[i].color() << "  " << pIn[i].anti_color();
  JSWARN << "pid = " << pIn[i].pid() << " E = " << pIn[i].e()
         << " px = " << pIn[i].p(1) << " py = " << pIn[i].p(2)
         << "  pz = " << pIn[i].p(3) << " virtuality = " << pIn[i].t()
         << " form_time in fm = " << pIn[i].form_time()
         << " split time = " << pIn[i].form_time() + pIn[i].x_in().t();
}

void Matter::DoEnergyLoss(double deltaT, double time, double Q2,
                          vector<Parton> &pIn, vector<Parton> &pOut) {

  if (std::isnan(pIn[0].e()) || std::isnan(pIn[0].px()) ||
      std::isnan(pIn[0].py()) || std::isnan(pIn[0].pz()) ||
      std::isnan(pIn[0].t()) || std::isnan(pIn[0].form_time())) {
    JSINFO << BOLDYELLOW << "Parton on entry busted on time step " << time;
    Dump_pIn_info(0, pIn);
  }


  VERBOSE(8) << MAGENTA << "SentInPartons Signal received : " << deltaT << " "
             << Q2 << " " << pIn.size();

  VERBOSE(8) << BOLDYELLOW
             << " ********************************************************** ";

  double rNum;

  double delT = deltaT;
  double Time = time * fmToGeVinv;
  double deltaTime = delT * fmToGeVinv;
  double ehat = 0;
  double ehat_over_T2 = 10.0;

  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  VERBOSE(8) << MAGENTA << " the time in fm is " << time
             << " The time in GeV-1 is " << Time;
  VERBOSE(8) << MAGENTA << "pid = " << pIn[0].pid() << " E = " << pIn[0].e()
             << " px = " << pIn[0].p(1) << " py = " << pIn[0].p(2)
             << "  pz = " << pIn[0].p(3) << " virtuality = " << pIn[0].t()
             << " form_time in fm = " << pIn[0].form_time()
             << " split time = " << pIn[0].form_time() + pIn[0].x_in().t();
  VERBOSE(8) << " color = " << pIn[0].color()
             << " anti-color = " << pIn[0].anti_color();

  unsigned int ShowerMaxColor = pIn[0].max_color();
  unsigned int CurrentMaxColor;

  if (pIn[0].max_color() < MaxColor) {
    pIn[0].set_max_color(MaxColor);

  } else {
    MaxColor = pIn[0].max_color();
  }

  //VERBOSE(8) << MAGENTA << " Max color = " << MaxColor;
  //JSDEBUG << " For Matter, the qhat in GeV^-3 = " << qhat ;

  double qhatbrick;
  if (brick_med)
    qhatbrick = qhat0 / 3.0;
  qhat = qhat0;

  VERBOSE(8) << " qhat0 = " << qhat0 << " qhat = " << qhat;

  for (int i = 0; i < pIn.size(); i++) {

    // Reject photons
    if (pIn[i].pid() == photonid) {
      if(pIn[i].pstat() != 22) {
            pIn[i].set_stat(22);
      	    VERBOSE(1) << BOLDYELLOW
                       << " A photon was RECEIVED with px = " << pIn[i].px()
                       << " from framework and sent back ";

            pOut.push_back(pIn[i]);
      }
      return;
    }

    // Reject photons

    if (std::abs(pIn[i].pstat()) == 1) {

      //          JSINFO << BOLDYELLOW << " A recoil was  RECEIVED with px = " << pIn[i].px() << " py = " << pIn[i].py() << " pz = " << pIn[i].pz() << " E = " << pIn[i].e() << " from framework and sent back " ;
      //          JSINFO << BOLDYELLOW << "t=" << " *  parton formation spacetime point= "<< pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z();
      //          Dump_pIn_info(i,pIn);

      //          pOut.push_back(pIn[i]);

      return;
    }

    VERBOSE(2) << BOLDYELLOW
               << " *  parton formation spacetime point= " << pIn[i].x_in().t()
               << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  "
               << pIn[i].x_in().z();

    //JSINFO << MAGENTA << " particle rest mass = " << pIn[i].restmass();

    int jet_stat = pIn[i].pstat(); // daughter of recoil will always be recoil

    JSDEBUG << MAGENTA << "Matter -- status: " << pIn[i].pstat()
            << "  energy: " << pIn[i].e() << " color: " << pIn[i].color()
            << "  " << pIn[i].anti_color() << "  clock: " << time;

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
    VERBOSE(2) << BOLDYELLOW << " velocityMod = " << velocityMod;

    if (pIn[i].form_time() < 0.0)
      pIn[i].set_jet_v(velocity); // jet velocity is set only once
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

    double SpatialRapidity =
        0.5 * std::log((now_R0 + now_Rz) / (now_R0 - now_Rz));

    //JSINFO << MAGENTA <<  " Particle Rapidity = " << SpatialRapidity ;

    if (std::isnan(velocityMod) || std::isnan(velocity[1]) ||
        std::isnan(velocity[2]) || std::isnan(velocity[3]) ||
        std::isinf(velocityMod) || std::isinf(velocity[1]) ||
        std::isinf(velocity[2]) || std::isinf(velocity[3])) {
      JSINFO << BOLDYELLOW << "Zeroth instance";
      JSINFO << BOLDYELLOW << "time, initR0, initRx, initRy, initRz=" << time
             << ", " << initR0 << ", " << initRx << ", " << initRy << ", "
             << initRz;
      JSINFO << BOLDYELLOW << "Vx, Vy, Vz =" << velocity[1] << ", "
             << velocity[2] << ", " << velocity[3];
      JSINFO << BOLDYELLOW << "initVx, initVy, initVz =" << initVx << ", "
             << initVy << ", " << initVz;
      Dump_pIn_info(i, pIn);
    }

    initEner = pIn[i].e(); // initial Energy of parton
    if (!in_vac) {
      if (GetJetSignalConnected())
        length = fillQhatTab(SpatialRapidity);
      else {
        JSWARN << "Couldn't find a hydro module attached!";
        throw std::runtime_error(
            "Please attach a hydro module or set in_vac to 1 in the XML file");
      }
    }
    if (brick_med)
      length = brick_length *
               fmToGeVinv; /// length in GeV-1 will have to changed for hydro
    //if(brick_med) length = 5.0*fmToGeVinv; /// length in GeV-1 will have to changed for hydro

    // SC
    zeta = ((xStart[0] + initRdotV) / std::sqrt(2)) * fmToGeVinv;

    VERBOSE(8) << BOLDYELLOW << " zeta = " << zeta;

    // if(now_R0^2-now_Ri^2<0) print out pIn info and exit

    if (std::isinf(now_R0) || std::isnan(now_R0) || std::isinf(now_Rz) ||
        std::isnan(now_Rz) || std::abs(now_Rz) > now_R0) {
      JSINFO << BOLDYELLOW << "First instance";
      JSINFO << BOLDYELLOW << "now_R for vector is:" << now_R0 << ", " << now_Rx
             << ", " << now_Ry << ", " << now_Rz;
      JSINFO << BOLDYELLOW << "time, initR0, initRx, initRy, initRz=" << time
             << ", " << initR0 << ", " << initRx << ", " << initRy << ", "
             << initRz;
      JSINFO << BOLDYELLOW << "initVx, initVy, initVz =" << initVx << ", "
             << initVy << ", " << initVz;
      JSINFO << BOLDYELLOW << "velocityMod=" << velocityMod;
      JSINFO << BOLDYELLOW << "initVMod="
             << std::sqrt(initVx * initVx + initVy * initVy + initVz * initVz);
      Dump_pIn_info(i, pIn);
      //exit(0);
    }
    double boostedTStart = tStart * cosh(SpatialRapidity);
    if (!in_vac && now_R0 >= boostedTStart) {
      if (now_R0 * now_R0 < now_Rz * now_Rz)
        cout << "Warning 1: " << now_R0 << "  " << now_Rz << endl;
      GetHydroCellSignal(now_R0, now_Rx, now_Ry, now_Rz, check_fluid_info_ptr);
      //VERBOSE(8)<<MAGENTA<<"Temperature from medium = "<<check_fluid_info_ptr->temperature;
      now_temp = check_fluid_info_ptr->temperature;
      //JSINFO << BOLDYELLOW << "Matter time = " << now_R0 << " x = " << now_Rx << " y = " << now_Ry << " z = " << now_Rz << " temp = " << now_temp;
      //JSINFO << BOLDYELLOW << "Matter initVx, initVy, initVz =" << initVx << ", " << initVy << ", " << initVz;
      //JSINFO << BOLDYELLOW << "Matter velocityMod=" << velocityMod;
    } else {
      now_temp = 0.0;
    }

    int pid = pIn[i].pid();

    if (pIn[i].form_time() <
        0.0) /// A parton without a virtuality or formation time, must set...
    {

      VERBOSE(8) << BOLDYELLOW << " pid = " << pIn[i].pid()
                 << " E = " << pIn[i].e();

      if ((pIn[i].t() < 0.0) &&
          ((pIn[i].form_time() < -0.1 - rounding_error) ||
           (pIn[i].form_time() > -0.1 + rounding_error))) {
        JSWARN << " parton with a negative virtuality was sent to Matter and "
                  "will now have its virtuality reset!, press 1 and return to "
                  "proceed... ";
        // cin >> blurb; //remove the input to prevent an error caused by heavy quark from pythia (by Chathuranga)
      }

      iSplit = 0; // (anti)quark splitting into (anti)quark + gluon
      if (pIn[i].pid() == gid) {
        JSDEBUG << " parton is a gluon ";
        iSplit = 1; // gluon splitting to two gluons
      } else {
        JSDEBUG << " parton is a quark ";
      }

      //tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, pIn[i].e()*pIn[i].e() ,zeta , iSplit);

      // SC:
      double pT2 = pIn[i].p(1) * pIn[i].p(1) + pIn[i].p(2) * pIn[i].p(2);
      double max_vir;
      if (vir_factor < 0.0)
        max_vir =
            std::abs(vir_factor) * (pIn[i].e() * pIn[i].e() - pIn[i].restmass() * pIn[i].restmass());
      else
        max_vir = pT2 * vir_factor;

      if (max_vir <= QS * QS) {
        tQ2 = 0.0;
      } else {
        JSDEBUG << BOLDYELLOW << " at x,y,z,t = " << pIn[i].x_in().x() << "  "
                << pIn[i].x_in().y() << "  " << pIn[i].x_in().z() << "  "
                << pIn[i].x_in().t();
        if (abs(pIn[i].pid()) == 4 || abs(pIn[i].pid()) == 5) {

          double min_vir =
              (QS * QS / 2.0) *
              (1.0 + std::sqrt(1.0 + 4.0 * pIn[i].restmass() *
                                         pIn[i].restmass() / QS / QS));
          if (max_vir > min_vir) {
            tQ2 =
                generate_vac_t_w_M(pIn[i].pid(), pIn[i].restmass(), pIn[i].nu(),
                                   QS * QS / 2.0, max_vir, zeta, iSplit);
          } else {
            tQ2 = QS * QS;
          }
          //  std::ofstream tdist;
          //  tdist.open("tdist_heavy.dat", std::ios::app);
          //  tdist << tQ2 << endl;
          //  tdist.close();

          VERBOSE(8) << BOLDYELLOW << " virtuality calculated as = " << tQ2;
        } else if (pIn[i].pid() == gid) {
          tQ2 = generate_vac_t_w_M(pIn[i].pid(), pIn[i].restmass(), pIn[i].nu(),
                                   QS * QS / 2.0, max_vir, zeta, iSplit);
        } else {
          tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS * QS / 2.0,
                               max_vir, zeta, iSplit);
        }
      }

      // SC: if Matter_on = false, set zero virtuality and Matter will not do parton shower
      if (matter_on) {
        pIn[i].set_t(tQ2); // Also resets momentum!
        VERBOSE(8) << BOLDYELLOW << " virtuality set to " << tQ2
                   << " max virtuality allowed = " << max_vir;
        VERBOSE(8) << BOLDYELLOW << " ID = " << pIn[i].pid()
                   << " Color = " << pIn[i].color()
                   << " Anti-Color = " << pIn[i].anti_color();
        VERBOSE(8) << BOLDYELLOW << " E = " << pIn[i].e()
                   << " px = " << pIn[i].px() << " py = " << pIn[i].py()
                   << " pz = " << pIn[i].pz();

      } else
        pIn[i].set_t(rounding_error);
      //else pIn[i].set_t(0.0);

      pIn[i].set_mean_form_time();
      double ft = generate_L(pIn[i].mean_form_time());
      pIn[i].set_form_time(ft);


      pIn[i].set_min_color(pIn[i].color());
      pIn[i].set_min_anti_color(pIn[i].anti_color());
      MaxColor = pIn[i].max_color();


      // VERBOSE OUTPUT ON INITIAL STATUS OF PARTICLE:
      VERBOSE(8);
      VERBOSE(8) << " *********************************************************"
                    "******************** ";
      VERBOSE(8) << BOLDYELLOW << " ID = " << pIn[i].pid()
                 << " Color = " << pIn[i].color()
                 << " Anti-Color = " << pIn[i].anti_color();
      VERBOSE(8) << BOLDYELLOW << " E = " << pIn[i].e()
                 << " px = " << pIn[i].px() << " py = " << pIn[i].py()
                 << " pz = " << pIn[i].pz();
      VERBOSE(8) << BOLDYELLOW << " *  New generated virtuality = " << tQ2
                 << " Mean formation time = " << pIn[i].mean_form_time();
      VERBOSE(8) << BOLDYELLOW << " *  set new formation time to "
                 << pIn[i].form_time();
      VERBOSE(8) << BOLDYELLOW << " * Maximum allowed virtuality = "
                 << pIn[i].e() * pIn[i].e() -
                        pIn[i].restmass() * pIn[i].restmass()
                 << "   Minimum Virtuality = " << QS * QS;
      VERBOSE(8) << " * Qhat = " << qhat << "  Length in fm = " << length / 5.0;
      VERBOSE(8) << " * Jet velocity = " << pIn[i].jet_v().comp(0) << " "
                 << pIn[i].jet_v().comp(1) << "  " << pIn[i].jet_v().comp(2)
                 << "  " << pIn[i].jet_v().comp(3);
      VERBOSE(8) << " * reset location of parton formation = "
                 << pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  "
                 << pIn[i].x_in().y() << "  " << pIn[i].x_in().z();
      VERBOSE(8) << " *********************************************************"
                    "******************** ";
      VERBOSE(8);
      // end VERBOSE OUTPUT:
    }

    // SC: Q0 can be changed based on different setups
    if (in_vac) { // for vaccuum
      qhat = 0.0;
      if (Q00 < 0.0)
        Q0 = 1.0; // set Q0 = 1 if Q00 < 0
      else
        Q0 = Q00;
    } else { // for medium
      double tempEner = initEner;
      qhat = fncQhat(zeta);
      ehat = 0.0;

      if (now_temp > 0.0)
        ehat = 0.0 * qhat / 4.0 / now_temp;
      VERBOSE(8) << BOLDYELLOW << "at Origin of parton, qhat = " << qhat
                 << " ehat = " << ehat;

      // set the minimum Q0 for Matter using Q^2 = qhat * tau  ELSE use the positive value in XML
      if (Q00 < 0.0) { // use dynamical Q0 if Q00 < 0
        if (pid == gid)
          Q0 = sqrt(sqrt(2.0 * tempEner * qhat * sqrt(2.0)));
        else
          Q0 = sqrt(sqrt(2.0 * tempEner * qhat * sqrt(2.0) / Ca * Cf));
        if (Q0 < 1.0)
          Q0 = 1.0;
        if (zeta > length)
          Q0 = 1.0;
      } else {
        Q0 = Q00;
      }
    }

    // I dont care what you say, we are not doing pQCD below 1 GeV
    if (Q0 < 1.0)
      Q0 = 1.0;

    //if (pIn[i].t() > QS + rounding_error)
    if (pIn[i].t() > Q0 * Q0 + rounding_error ||
        ((!in_vac) && now_temp <= T0 &&
         pIn[i].t() > QS * QS + rounding_error)) {

      TakeResponsibilityFor(
          pIn[i]); // Generate error if another module already has responsibility.
      double decayTime = pIn[i].mean_form_time();

      //JSDEBUG << "  deltaT = " << deltaT;
      VERBOSE(8) << " parton origin time = " << pIn[i].x_in().t()
                 << " parton formation time = " << pIn[i].form_time();
      VERBOSE(8) << " parton id " << pIn[i].pid()
                 << " parton virtuality = " << pIn[i].t();
      //JSDEBUG << " parton momentum " << pIn[i].e() << "  " << pIn[i].px() << "  " << pIn[i].py() << "  " << pIn[i].pz();

      double splitTime = pIn[i].form_time() + pIn[i].x_in().t();

      VERBOSE(8) << " splitTime = " << splitTime;
      VERBOSE(8) << " qhat before splitime loop = " << qhat;

      if (splitTime <
          time) // it is time to split and calculate the effect of scattering
      {

        VERBOSE(8) << "SPLIT in Matter";

        // SC: add elastic scattering that generates recoiled and back-reaction partons
        double el_dt = 0.1;
        double el_p0[5];
        double el_CR;
        double el_rand;
        double HQ_mass;

        if (pid == gid)
          el_CR = Ca;
        else
          el_CR = Cf;

        for (int j = 1; j <= 3; j++)
          el_p0[j] = pIn[i].p(j);
        el_p0[0] = pIn[i].e();
        el_p0[4] = pIn[i].t();
        //HQ_mass = 1.27;////1.5 given by pythia originally , 1.27 by Wikipedia
        HQ_mass = pIn[i].restmass();

        for (double el_time = initR0; el_time < time + rounding_error;
             el_time = el_time + el_dt) {

          if (in_vac)
            continue;
          if (!recoil_on)
            continue;

          double boostedTStart = tStart * std::cosh(SpatialRapidity);
          if (el_time < boostedTStart)
            continue;

          if (abs(pid) == 5)
            continue; // recoil not ready for b quark yet

          double el_rx = initRx + (el_time - initR0) * initVx;
          double el_ry = initRy + (el_time - initR0) * initVy;
          double el_rz = initRz + (el_time - initR0) * initVz;

          double tempLoc, sdLoc, vxLoc, vyLoc, vzLoc, qhatLoc, enerLoc;
          double betaLoc, gammaLoc, flowFactor;
          int hydro_ctl;
          double dt_lrf;

          double pc0[4] = {0.0};
          double vc0[4] = {0.0};
          double soln_alphas, prob_el;
          double muD2;
          double recordE0;

          // Convert hard parton momentum to onshell
          pc0[1] = el_p0[1];
          pc0[2] = el_p0[2];
          pc0[3] = el_p0[3];
          if (abs(pid) == 4 || abs(pid) == 5)

            pc0[0] = sqrt(pc0[1] * pc0[1] + pc0[2] * pc0[2] + pc0[3] * pc0[3] + 1.27*1.27);
                         // HQ_mass * HQ_mass);
          else
            pc0[0] = sqrt(pc0[1] * pc0[1] + pc0[2] * pc0[2] + pc0[3] * pc0[3]);

          recordE0 = pc0[0];

          if (std::isinf(el_time) || std::isnan(el_time) || std::isinf(el_rz) ||
              std::isnan(el_rz) || std::abs(el_rz) > el_time) {
            JSWARN << "Second instance";
            JSWARN << "el_vector for vector is:" << el_time << ", " << el_rx
                   << ", " << el_ry << ", " << el_rz;
            JSWARN << "initR0, initRx, initRy, init Rz=" << initR0 << ", "
                   << initRx << ", " << initRy << ", " << initRz;
            JSWARN << "initVx, initVy, initVz =" << initVx << ", " << initVy
                   << ", " << initVz;
            JSWARN << "velocityMod=" << std::setprecision(20) << velocityMod;
            JSWARN << "initVMod=" << std::setprecision(20)
                   << std::sqrt(initVx * initVx + initVy * initVy +
                                initVz * initVz);
            Dump_pIn_info(i, pIn);
            //exit(0);
          }
          GetHydroCellSignal(el_time, el_rx, el_ry, el_rz,
                             check_fluid_info_ptr);
          VERBOSE(8) << MAGENTA << "Temperature from medium = "
                     << check_fluid_info_ptr->temperature;

          tempLoc = check_fluid_info_ptr->temperature;
          sdLoc = check_fluid_info_ptr->entropy_density;
          vxLoc = check_fluid_info_ptr->vx;
          vyLoc = check_fluid_info_ptr->vy;
          vzLoc = check_fluid_info_ptr->vz;

          vc0[1] = vxLoc;
          vc0[2] = vyLoc;
          vc0[3] = vzLoc;

          hydro_ctl = 0;
          if (hydro_ctl == 0 && tempLoc >= hydro_Tc) {

            trans(vc0, pc0);
            enerLoc = pc0[0];
            transback(vc0, pc0);

            betaLoc = sqrt(vxLoc * vxLoc + vyLoc * vyLoc + vzLoc * vzLoc);
            gammaLoc = 1.0 / sqrt(1.0 - betaLoc * betaLoc);
            flowFactor =
                gammaLoc *
                (1.0 - (initVx * vxLoc + initVy * vyLoc + initVz * vzLoc));

	    //if(run_alphas==1){ alphas= 4*pi/(9.0*log(2*enerLoc*tempLoc/0.04));}

	    /*            if (qhat0 < 0.0) { // calculate qhat with alphas
              muD2 = 6.0 * pi * alphas * tempLoc * tempLoc;
              if (enerLoc > 2.0 * pi * tempLoc)
                qhatLoc = Ca * 50.4864 / pi * pow(alphas, 2) * pow(tempLoc, 3) *
                          log(5.7 * enerLoc * tempLoc / 4.0 / muD2);
              else
                qhatLoc = Ca * 50.4864 / pi * pow(alphas, 2) * pow(tempLoc, 3) *
                          log(5.7 * 2.0 * pi * tempLoc * tempLoc / 4.0 / muD2);
              if (qhatLoc < 0.0)
                qhatLoc = 0.0;
            } else { // use input qhat
              if (brick_med)
                qhatLoc = qhat0 * 0.1973;
              else
                qhatLoc = qhat0 / 96.0 * sdLoc * 0.1973;
		}*/

	    //GeneralQhatFunction(int QhatParametrizationType, double Temperature, double EntropyDensity, double FixAlphas,  double Qhat0, double E, double muSquare)
	    double muSquare= pIn[i].t(); //Virtuality of the parent; Revist this when q-hat is virtuality dependent
	    qhatLoc= GeneralQhatFunction(QhatParametrizationType, tempLoc, sdLoc, alphas, qhat0, enerLoc, muSquare);
          } else { // outside the QGP medium
            continue;
          }

          // Calculating the delta t in the fluid rest frame
          if (el_time + el_dt <= time)
            dt_lrf = el_dt * flowFactor;
          else
            dt_lrf = (time - el_time) * flowFactor;

	  //if(run_alphas==1){ alphas= 4*pi/(9.0*log(2*enerLoc*tempLoc/0.04));}

          // solve alphas
          if (qhat0 < 0.0 || QhatParametrizationType==0 || QhatParametrizationType==1 || QhatParametrizationType==5 ||
              QhatParametrizationType==6 || QhatParametrizationType==7 )
            soln_alphas = alphas;
          else
            soln_alphas = solve_alphas(qhatLoc, enerLoc, tempLoc);

          // Calculate the proability of elastic scattering in time delta t in fluid rest frame
          muD2 = 6.0 * pi * soln_alphas * tempLoc * tempLoc;
          prob_el = 42.0 * zeta3 * el_CR  * tempLoc / 6.0 / pi /
                    pi * dt_lrf / 0.1973;

	        prob_el=prob_el*ModifiedProbability(QhatParametrizationType, tempLoc, sdLoc, enerLoc, pIn[i].t());

          el_rand = ZeroOneDistribution(*GetMt19937Generator());

          //cout << "  qhat: " << qhatLoc << "  alphas: " << soln_alphas << "  ener: " << enerLoc << "  prob_el: " << prob_el << "  " << el_rand << endl;

          if (el_rand < prob_el) { // elastic scattering happens
            //cout << "elastic scattering happens" << endl;
            int CT = -1; //could be removed
            int pid0 = pid;  //ivan
            int pid2 = -999;  //ivan
            int pid3 = -999; //ivan

            double pc2[4] = {0.0}; // final recoil thermal parton //can be removed
            double pc3[4] = {0.0}; // initial thermal parton  //can be removed
            double pc4[4] = {0.0}; // not used //can be removed
            double qt = 0.0;       // not used //not used
            unsigned int el_max_color, el_color0, el_anti_color0, el_color2,
                el_anti_color2, el_color3, el_anti_color3; //not used

            el_max_color = pIn[i].max_color(); //not used
            el_color0 = pIn[i].color(); //not used
            el_anti_color0 = pIn[i].anti_color(); //not used

            //pid0 = pid; //pid0 no use

            // deterimine channel
            //flavor(CT,pid0,pid2,pid3);
            trans(vc0,pc0); //ivan
            elasticCollision.elastic_kinematics(tempLoc,pid0,pid2,pid3,pc0,pc2,pc3); //ivan
            transback(vc0,pc0);   //leading ivan
            transback(vc0,pc2);   //recoil ivan
            transback(vc0,pc3);   //hole ivan
            
            //flavor(CT, pid0, pid2, pid3, el_max_color, el_color0,
            //       el_anti_color0, el_color2, el_anti_color2, el_color3,
            //       el_anti_color3);

            //cout << "color: " << el_color0 << "  " << el_anti_color0 << "  " << el_color2 << "  " << el_anti_color2 << "  " << el_color3 << "  " << el_anti_color3 << "  max: " << el_max_color << endl;

            // do scattering
            //if (CT == 11 || CT == 12) { // for heavy quark scattering
            //  collHQ22(CT, tempLoc, muD2, vc0, pc0, pc2, pc3, pc4, qt);
            //} else { // for light parton scattering
            //  colljet22(CT, tempLoc, muD2, vc0, pc0, pc2, pc3, pc4, qt);
            //}

            if (pc0[0] < pc2[0] && abs(pid0) != 4 &&
                pid0 ==
                    pid2) { //disable switch for heavy quark, only allow switch for identical particles
              double p0temp[4] = {0.0};
              for (int k = 0; k <= 3; k++) {
                p0temp[k] = pc2[k];
                pc2[k] = pc0[k];
                pc0[k] = p0temp[k];
              }
            }

            //Amit: do not add recoil parton with E=0,Px=0,Py=0,Pz=0 into pOut vector
            if (pc2[0] < rounding_error || pc3[0] < rounding_error)
              continue;

            // push out recoied and back-reaction (negative) partons
            // need to add color information later!!!
            double el_vertex[4];
            el_vertex[0] = el_time;
            el_vertex[1] = el_rx;
            el_vertex[2] = el_ry;
            el_vertex[3] = el_rz;

            int iout;
            double ft;

            pc2[0] = sqrt(pc2[1] * pc2[1] + pc2[2] * pc2[2] + pc2[3] * pc2[3] +
                          rounding_error);

            if (std::isnan(pc2[1]) || std::isnan(pc2[2]) ||
                std::isnan(pc2[3]) || std::isinf(pc2[1]) ||
                std::isinf(pc2[2]) || std::isinf(pc2[3])) {
              JSWARN << "recoil in Matter instance 1: pc[0]=" << pc2[0]
                     << ", pc2[1]=" << pc2[1] << ", pc2[2]=" << pc2[2]
                     << ", pc2[3]=" << pc2[3];
            }

            pOut.push_back(Parton(0, pid2, 1, pc2, el_vertex)); // recoiled
            iout = pOut.size() - 1;
            pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
            pOut[iout].set_mean_form_time();
            ft = 10000.0; /// a really large formation time.
            pOut[iout].set_form_time(ft);
            ////pOut[iout].set_color(el_color2);
            ////pOut[iout].set_anti_color(el_anti_color2);
            ////pOut[iout].set_max_color(max_color);
            ////pOut[iout].set_min_color(pIn[i].min_color());
            ////pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
            ////// comment out realistic color above, assume colorless for recoiled and back-reaction parton
            ////// for the convenience of color string fragmentation in Pythia
            pOut[iout].set_color(0);
            pOut[iout].set_anti_color(0);
            pOut[iout].set_max_color(pIn[i].max_color());
            pOut[iout].set_min_color(pIn[i].min_color());
            pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

            pOut.push_back(
                Parton(0, pid3, -1, pc3, el_vertex)); // back reaction
            iout = pOut.size() - 1;
            pOut[iout].set_jet_v(velocity_jet);
            pOut[iout].set_mean_form_time();
            ft = 10000.0; /// STILL a really large formation time.
            pOut[iout].set_form_time(ft);
            ////pOut[iout].set_color(el_color3);
            ////pOut[iout].set_anti_color(el_anti_color3);
            ////pOut[iout].set_max_color(max_color);
            ////pOut[iout].set_min_color(pIn[i].min_color());
            ////pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
            pOut[iout].set_color(0);
            pOut[iout].set_anti_color(0);
            pOut[iout].set_max_color(pIn[i].max_color());
            pOut[iout].set_min_color(pIn[i].min_color());
            pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

            ////// assumption: pIn doesn't change color in elastic scattering
            ////pIn[i].set_color(el_color0);
            ////pIn[i].set_anti_color(el_anti_color0);
            ////pIn[i].set_max_color(el_max_color);

          } // end if scattering

          // E1'+E2 = E3'+E4 (all massless in scattering)
          // Now I want E1+E2 = E3+E4 where 1 and 3 have mass
          // -> E1-E1' = E3-E3' -> E3 = E1-E1'+E3'
          // m3^2 = E3^2-E3'^2
          for (int j = 1; j <= 3; j++)
            el_p0[j] = pc0[j];

          el_p0[0] = el_p0[0] - recordE0 + pc0[0];
          el_p0[4] = el_p0[0] * el_p0[0] - pc0[0] * pc0[0];
          if (el_p0[4] < 0)
            cout << "complain negative virt" << endl;

        } // end time loop for elastic scattering

        // do split
        double t_used = pIn[i].t();
        //if (t_used<QS)  t_used = QS; // SC: not necessary
        double tau_form = 2.0 * pIn[i].nu() / t_used;
        double z_low = QS * QS / t_used / 2.0;
        double z_hi = 1.0 - z_low;

        VERBOSE(8) << " zeta = " << zeta;

        if (pid == gid) { // gluon
          double val1 =
              SplittingFunction::P_z_gg_int(z_low, z_hi, zeta, t_used, tau_form, pIn[i].nu(),*this);
              //JSINFO<<" "<<val1;
          double val2 =
              nf * SplittingFunction::P_z_qq_int(z_low, z_hi, zeta, t_used, tau_form, pIn[i].nu(),*this);
          double M = PythiaFunction.particleData.m0(cid);
          z_low = (QS * QS + 2.0 * M * M) / t_used / 2.0;
          z_hi = 1.0 - z_low;
          double val3 = 0.0;
          if (z_hi > z_low)
            val3 = SplittingFunction::P_z_qq_int_w_M_vac_only(M, z_low, z_hi, zeta, t_used,
                                           tau_form, pIn[i].nu(),*this);
          if (t_used < (QS * QS + 2.0 * M * M))
            val3 = 0.0;

          M = PythiaFunction.particleData.m0(bid);
          z_low = (QS * QS + 2.0 * M * M) / t_used / 2.0;
          z_hi = 1.0 - z_low;
          double val4 = 0.0;
          if (z_hi > z_low)
            val4 = SplittingFunction::P_z_qq_int_w_M_vac_only(M, z_low, z_hi, zeta, t_used,
                                           tau_form, pIn[i].nu(),*this);
          if (t_used < (QS * QS + 2.0 * M * M))
            val4 = 0.0;

          VERBOSE(8) << BOLDYELLOW << " val1 = " << val1 << " val2 = " << val2
                     << " val3 = " << val3 << " val4 = " << val4;

          if (val1 < 0.0 || val2 < 0.0 || val3 < 0.0 || val4 < 0.0) {
            cerr << " minus log of sudakov negative val1 , val2 , val3, val4 = "
                 << val1 << "  " << val2 << "  " << val3 << "  " << val4
                 << endl;
            throw std::runtime_error("minus log of sudakov negative");
            // cin >> blurb ;
          }

          double ratio1 = val1 / (val1 + val2 + val3 + val4);
          double ratio2 = (val1 + val2) / (val1 + val2 + val3 + val4);
          double ratio3 = (val1 + val2 + val3) / (val1 + val2 + val3 + val4);
          double r = ZeroOneDistribution(*GetMt19937Generator());
          if (r >= ratio1 && r < ratio2) { // qqbar

            double r2 = ZeroOneDistribution(*GetMt19937Generator());

            // assign flavors
            if (r2 > 0.6666) {
              pid_a = uid;
              pid_b = -1 * uid;
            } else if (r2 > 0.3333) {
              pid_a = did;
              pid_b = -1 * did;

            } else {
              pid_a = sid;
              pid_b = -1 * sid;
            }
            iSplit = 2;
          } else if (r >= ratio2 && r < ratio3) {
            pid_a = cid;
            pid_b = -cid;
            iSplit = 4;

            VERBOSE(1) << BOLDYELLOW << " split to c c-bar";

          } else if (r >= ratio3) {
            pid_a = bid;
            pid_b = -bid;
            iSplit = 5;

            VERBOSE(1) << BOLDYELLOW << " Split to b b-bar ";
          } else { // gg
            pid_a = gid;
            pid_b = gid;
            iSplit = 1;
          }
        } else if (std::abs(pIn[i].pid()) < 4) { // we had a light quark

          double relative_charge = 0.0;
          if ((std::abs(pIn[i].pid()) > 0) && (std::abs(pIn[i].pid()) < 4))
            relative_charge = 1.0 / 9.0;
          if (std::abs(pIn[i].pid()) == 2)
            relative_charge = relative_charge * 4.0;
          if (std::abs(pIn[i].pid()) == 4)
            relative_charge = 4.0 / 9.0;
          if (std::abs(pIn[i].pid()) == 5)
            relative_charge = 1.0 / 9.0;

          double ProbGluon =
              1.0 - SudakovFunction::sudakov_Pqg(QS * QS / 2, pIn[i].t(), zeta, pIn[i].nu(),*this);
          double ProbPhoton =
              1.0 -
              std::pow(SudakovFunction::sudakov_Pqp(QS * QS / 2, pIn[i].t(), zeta, pIn[i].nu(),*this),
                       relative_charge);

          double val = ProbGluon / (ProbGluon + ProbPhoton);

          VERBOSE(8) << MAGENTA
                     << " probability of gluon radiation from quark = " << val;

          double r2 = ZeroOneDistribution(*GetMt19937Generator());

          if (r2 <= val) { // light quark decay to quark and gluon
            pid_a = pid;
            pid_b = gid;
            iSplit = 0; // iSplit for quark radiating a gluon
          } else {      // light quark decay to quark and photon
            pid_a = pid;
            pid_b = photonid;
            iSplit = 3; // iSplit for quark radiating a photon
            photon_brem = true;
          }

        } else {
          // we had a heavy quark
          pid_a = pid;
          pid_b = gid;
          iSplit = 0;
        }

        // daughter virtualities
        double tQd1 = QS * QS;
        double tQd2 = QS * QS;

        double new_parent_p0 = el_p0[0];
        double new_parent_px = el_p0[1];
        double new_parent_py = el_p0[2];
        double new_parent_pz = el_p0[3];
        double new_parent_t = el_p0[4];
        double new_parent_pl = (new_parent_px * pIn[i].jet_v().x() +
                                new_parent_py * pIn[i].jet_v().y() +
                                new_parent_pz * pIn[i].jet_v().z()) /
                               mod_jet_v;
        if (new_parent_pl < 0.0) {
          VERBOSE(8) << BOLDYELLOW
                     << " parton traversing opposite to jet direction ... Just "
                        "letting you know ! ";
          // cin >> blurb ;
        }
        //JSINFO << BOLDYELLOW << " old virtuality = " << pIn[i].t() << " new virtuality = " << new_parent_t ;
        //double new_parent_pl = sqrt(pow(new_parent_px,2)+pow(new_parent_py,2)+pow(new_parent_pz,2));
        //double new_parent_vx = new_parent_px/new_parent_pl;
        //double new_parent_vy = new_parent_py/new_parent_pl;
        //double new_parent_vz = new_parent_pz/new_parent_pl;
        double new_parent_nu = (new_parent_p0 + new_parent_pl) / sqrt(2.0);

        //set color of daughters here
        unsigned int d1_col, d1_acol, d2_col, d2_acol, color, anti_color;
        //std::uniform_int_distribution<short> uni(101,103);
        //color = pIn[i].color();

        if (iSplit != 3) // not photon radiation, generate new colors
        {
          max_color = pIn[i].max_color();
          //if (pIn[i].anti_color()>maxcolor) color = pIn[i].anti_color();
          JSDEBUG << " old max color = " << max_color;
          max_color = ++MaxColor;
          color = max_color;
          anti_color = max_color;
          pIn[i].set_max_color(max_color);
          JSDEBUG << " new color = " << color;
        }

        if (iSplit == 1) ///< gluon splits into two gluons
        {
          d1_col = pIn[i].color();
          d2_col = color;
          d1_acol = anti_color;
          d2_acol = pIn[i].anti_color();
        } else if (
            iSplit ==
            0) ///< (anti-)quark splits into (anti-)quark + gluon, covers both light and heavy quarks (anti-quarks)
        {
          if (pIn[i].pid() > 0) // parent is a quark
          {
            d1_col = color;
            d1_acol = 0;
            d2_col = pIn[i].color();
            d2_acol = anti_color;
          } else {
            d1_col = 0; // parent is an anti-quark
            d1_acol = anti_color;
            d2_col = color;
            d2_acol = pIn[i].anti_color();
          }
        } else if (
            iSplit == 2 || iSplit == 4 ||
            iSplit ==
                5) // gluon splits into quark anti-quark, c anti-c, b anti-b
        {
          d1_col = pIn[i].color();
          d1_acol = 0;
          d2_acol = pIn[i].anti_color();
          d2_col = 0;
        } else if (
            iSplit ==
            3) // radiating a photon has col = acol = 0, all color remains in quark(anti-quark)
        {
          d1_col = pIn[i].color();
          d1_acol = pIn[i].anti_color();
          d2_col = 0;
          d2_acol = 0;
        } else {
          throw std::runtime_error("error in iSplit");
        }

        double l_perp2 = -1.0; // SC: initialization
        int ifcounter = 0;

        while ((l_perp2 <= Lambda_QCD * Lambda_QCD) && (ifcounter < 100)) {

          // if(abs(pid) == 4 || abs(pid) == 5)
          // {
          z = generate_vac_z_w_M(pid, pIn[i].restmass(), QS * QS / 2.0,
                                 pIn[i].t(), zeta, pIn[i].nu(), iSplit);
          // }
          // else if (pid == 21 && ( iSplit == 4 || iSplit == 5 ) )
          // {
          //     z = generate_vac_z(pid, QS/2.0, pIn[i].t(), zeta, pIn[i].nu(), iSplit);
          // }
          VERBOSE(8) << MAGENTA << " generated z = " << z;

          int iSplit_a = 0;
          if (pid_a == gid)
            iSplit_a = 1;

          // use pIn information to sample z above, but use new_parent to calculate daughter partons below

          if (std::abs(pid_a) == 4 || std::abs(pid_a) == 5) {
            double M = PythiaFunction.particleData.m0(pid_a);
            if (QS * QS * (1.0 + std::sqrt(1.0 + 4.0 * M * M / QS / QS)) / 2.0 <
                z * z * new_parent_t) {
              tQd1 = generate_vac_t_w_M(
                  pid_a, pIn[i].restmass(), z * new_parent_nu, QS * QS / 2.0,
                  z * z * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_a);

            } else {
              tQd1 = z * z * new_parent_t;
            }
          } else if (pid_a == 21) {
            double M = 0.0;
            if ((QS * QS) < z * z * new_parent_t) {
              tQd1 = generate_vac_t_w_M(
                  pid_a, M, z * new_parent_nu, QS * QS / 2.0,
                  z * z * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_a);
            } else {
              tQd1 = z * z * new_parent_t;
            }
            VERBOSE(8) << BOLDYELLOW << " daughter 1 virt generated = " << tQd1;

          } else // light quark anti-quark
          {
            if (z * z * new_parent_t > QS * QS) {
              tQd1 = generate_vac_t(
                  pid_a, z * new_parent_nu, QS * QS / 2.0, z * z * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_a);
            } else { // SC
              tQd1 = z * z * new_parent_t;
            }
          }

          int iSplit_b = 0;
          if (pid_b == gid)
            iSplit_b = 1;

          if (std::abs(pid_b) == 4 || std::abs(pid_b) == 5) {
            double M = PythiaFunction.particleData.m0(pid_b);

            if (QS * QS * (1.0 + std::sqrt(1.0 + 4.0 * M * M / QS / QS)) / 2.0 <
                (1.0 - z) * (1.0 - z) * new_parent_t) {
              tQd2 = generate_vac_t_w_M(
                  pid_b, pIn[i].restmass(), (1.0 - z) * new_parent_nu,
                  QS * QS / 2.0, (1.0 - z) * (1.0 - z) * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_b);

            } else {
              tQd2 = (1.0 - z) * (1.0 - z) * new_parent_t;
            }

          } else if (pid_b == 21) {
            double M = 0.0;
            if (QS * QS < (1.0 - z) * (1.0 - z) * new_parent_t) {
              tQd2 = generate_vac_t_w_M(
                  pid_b, M, (1.0 - z) * new_parent_nu, QS * QS / 2.0,
                  (1.0 - z) * (1.0 - z) * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_b);

            } else {
              tQd2 = (1.0 - z) * (1.0 - z) * new_parent_t;
            }
            VERBOSE(8) << BOLDYELLOW << " daughter virt generated = " << tQd2;

          } else {
            if (((1.0 - z) * (1.0 - z) * new_parent_t > QS * QS) &&
                (iSplit != 3)) {
              tQd2 = generate_vac_t(
                  pid_b, (1.0 - z) * new_parent_nu, QS * QS / 2.0,
                  (1.0 - z) * (1.0 - z) * new_parent_t,
                  zeta + std::sqrt(2) * pIn[i].form_time() * fmToGeVinv,
                  iSplit_b);
            } else { // SC
              tQd2 = (1.0 - z) * (1.0 - z) * new_parent_t;
            }

            if (iSplit == 3) {
              tQd2 = rounding_error; // forcing the photon to have no virtuality
            }
          }

          l_perp2 = new_parent_t * z * (1.0 - z) - tQd2 * z -
                    tQd1 * (1.0 - z); ///< the transverse momentum squared

          if (abs(pid) == 4 || abs(pid) == 5) {
            l_perp2 = new_parent_t * z * (1.0 - z) - tQd2 * z -
                      tQd1 * (1.0 - z) -
                      std::pow((1.0 - z) * pIn[i].restmass(),
                               2); ///< the transverse momentum squared
          } else if ((pid == 21) &&
                     (iSplit > 3)) // gluon decay into heavy quark anti-quark
          {

            double M = PythiaFunction.particleData.m0(pid_a);

            l_perp2 = new_parent_t * z * (1.0 - z) - tQd2 * z -
                      tQd1 * (1.0 - z) - M * M;
          } else {
            l_perp2 = new_parent_t * z * (1.0 - z) - tQd2 * z -
                      tQd1 * (1.0 - z); ///< the transverse momentum squared
          }

          ifcounter++;
        }

        // std::ofstream zdist;
        // zdist.open("zdist_heavy.dat", std::ios::app);
        // std::ofstream tdist;
        // tdist.open("tdist_heavy.dat", std::ios::app);
        // if (std::abs(pid_a) == 4 || std::abs(pid_a) == 5)
        // {
        //     tdist << tQd1<< endl;
        //     zdist << z << endl;
        // }
        // else if (std::abs(pid_b) == 4 || std::abs(pid_b) == 5)
        // {
        //     tdist << tQd2<< endl;
        //     zdist << z << endl;
        // }
        // tdist.close();
        // zdist.close();

        if (l_perp2 <= Lambda_QCD * Lambda_QCD)
          l_perp2 = Lambda_QCD * Lambda_QCD; ///< test if negative
        double l_perp = std::sqrt(
            l_perp2); ///< the momentum transverse to the parent parton direction
        VERBOSE(8) << BOLDYELLOW << " after ifcounter = " << ifcounter
                   << " l_perp2 = " << l_perp2;
        VERBOSE(8) << BOLDYELLOW << " z = " << z << " tQd1 = " << tQd1
                   << " tQd2 = " << tQd2;
        VERBOSE(8) << BOLDYELLOW << " pid_a = " << pid_a
                   << " pid_b = " << pid_b;

        // axis of split
        double angle = generate_angle();

        // KK: changed to x,y,z
        //double parent_perp = std::sqrt( pow(pIn[i].px(),2) + pow(pIn[i].py(),2) + pow(pIn[i].pz(),2) - pow(pIn[i].pl(),2) );
        //double mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
        double c_t =
            pIn[i].jet_v().z() /
            mod_jet_v; // c_t is cos(theta) for the jet_velocity unit vector
        double s_t = std::sqrt(1.0 - c_t * c_t); // s_t is sin(theta)

        double s_p = pIn[i].jet_v().y() / std::sqrt(pow(pIn[i].jet_v().x(), 2) +
                                                    pow(pIn[i].jet_v().y(), 2));
        double c_p = pIn[i].jet_v().x() / std::sqrt(pow(pIn[i].jet_v().x(), 2) +
                                                    pow(pIn[i].jet_v().y(), 2));

        VERBOSE(8) << BOLDYELLOW
                   << " Jet direction w.r.t. beam: theta = " << std::acos(c_t)
                   << " phi = " << std::acos(c_p);

        //double c_t = new_parent_vz;
        //double s_t = std::sqrt( 1.0 - c_t*c_t) ;
        //
        //double s_p = new_parent_vy/std::sqrt( pow(new_parent_vx,2) + pow(new_parent_vy,2) ) ;
        //double c_p = new_parent_vx/std::sqrt( pow(new_parent_vx,2) + pow(new_parent_vy,2) ) ;

        // First daughter
        double k_perp1[4];
        k_perp1[0] = 0.0;
        k_perp1[1] = z * (new_parent_px - new_parent_pl * s_t * c_p) +
                     l_perp * std::cos(angle) * c_t * c_p -
                     l_perp * std::sin(angle) * s_p;
        k_perp1[2] = z * (new_parent_py - new_parent_pl * s_t * s_p) +
                     l_perp * std::cos(angle) * c_t * s_p +
                     l_perp * std::sin(angle) * c_p;
        k_perp1[3] = z * (new_parent_pz - new_parent_pl * c_t) -
                     l_perp * std::cos(angle) * s_t;
        double k_perp1_2 =
            pow(k_perp1[1], 2) + pow(k_perp1[2], 2) + pow(k_perp1[3], 2);

        double M = 0.0;

        if ((std::abs(pid_a) == 4) || (std::abs(pid_a) == 5))
          M = PythiaFunction.particleData.m0(pid_a);

        double energy = (z * new_parent_nu + (tQd1 + k_perp1_2 + M * M) /
                                                 (2.0 * z * new_parent_nu)) /
                        std::sqrt(2.0);
        double plong = (z * new_parent_nu - (tQd1 + k_perp1_2 + M * M) /
                                                (2.0 * z * new_parent_nu)) /
                       std::sqrt(2.0);
        if (energy < 0.0) {
          JSWARN << " Energy negative after rotation, press 1 and return to "
                    "continue ";
          cin >> blurb;
        }

        double newp[4];
        newp[0] = energy;
        newp[1] = plong * s_t * c_p + k_perp1[1];
        newp[2] = plong * s_t * s_p + k_perp1[2];
        newp[3] = plong * c_t + k_perp1[3];

        VERBOSE(8) << MAGENTA << " D1 px = " << newp[1] << " py = " << newp[2]
                   << " pz = " << newp[3] << " E = " << newp[0];
        double newx[4];
        //newx[0] = time + deltaT;
        newx[0] = time;
        for (int j = 1; j <= 3; j++) {
          //newx[j] = pIn[i].x_in().comp(j) + (time + deltaT - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
          newx[j] = pIn[i].x_in().comp(j) +
                    (time - pIn[i].x_in().comp(0)) * velocity[j];
        }

        VERBOSE(8) << BOLDRED << " PiD - a = " << pid_a;

        pOut.push_back(Parton(0, pid_a, jet_stat, newp, newx));
        int iout = pOut.size() - 1;

        if (std::isnan(newp[1]) || std::isnan(newp[2]) || std::isnan(newp[3])) {
          JSINFO << MAGENTA << plong << " " << s_t << " " << c_p << " "
                 << k_perp1[1];

          JSINFO << MAGENTA << newp[0] << " " << newp[1] << " " << newp[2]
                 << " " << newp[3];
          cin >> blurb;
        }
        pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
        pOut[iout].set_mean_form_time();
        double ft = generate_L(pOut[iout].mean_form_time());
        pOut[iout].set_form_time(ft);
        pOut[iout].set_color(d1_col);
        pOut[iout].set_anti_color(d1_acol);
        pOut[iout].set_max_color(max_color);
        pOut[iout].set_min_color(pIn[i].min_color());
        pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

        VERBOSE(8) << BOLDRED << " virtuality of D 1 = " << pOut[iout].t();
        VERBOSE(8) << BOLDRED << " mass of parton = " << pOut[iout].restmass();

        // Second daughter
        double k_perp2[4];
        k_perp2[0] = 0.0;
        k_perp2[1] = (1.0 - z) * (new_parent_px - new_parent_pl * s_t * c_p) -
                     l_perp * std::cos(angle) * c_t * c_p +
                     l_perp * std::sin(angle) * s_p;
        k_perp2[2] = (1.0 - z) * (new_parent_py - new_parent_pl * s_t * s_p) -
                     l_perp * std::cos(angle) * c_t * s_p -
                     l_perp * std::sin(angle) * c_p;
        k_perp2[3] = (1.0 - z) * (new_parent_pz - new_parent_pl * c_t) +
                     l_perp * std::cos(angle) * s_t;
        double k_perp2_2 =
            pow(k_perp2[1], 2) + pow(k_perp2[2], 2) + pow(k_perp2[3], 2);

        M = 0.0;
        if ((std::abs(pid_b) == 4) || (std::abs(pid_b) == 5))
          M = PythiaFunction.particleData.m0(pid_b);
        ;

        energy =
            ((1.0 - z) * new_parent_nu +
             (tQd2 + k_perp2_2 + M * M) / (2.0 * (1.0 - z) * new_parent_nu)) /
            std::sqrt(2.0);
        plong =
            ((1.0 - z) * new_parent_nu -
             (tQd2 + k_perp2_2 + M * M) / (2.0 * (1.0 - z) * new_parent_nu)) /
            std::sqrt(2.0);

        //parent_perp = std::sqrt( pow(pIn[i].p(1),2) + pow(pIn[i].p(2),2) + pow(pIn[i].p(3),2) - pow(pIn[i].pl(),2) );
        //mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;

        if (energy < 0.0) {
          JSWARN << " Energy of 2nd daughter negative after rotation, press 1 "
                    "and return to continue ";
          cin >> blurb;
        }

        newp[0] = energy;
        newp[1] = plong * s_t * c_p + k_perp2[1];
        newp[2] = plong * s_t * s_p + k_perp2[2];
        newp[3] = plong * c_t + k_perp2[3];

        if (std::isnan(newp[1]) || std::isnan(newp[2]) || std::isnan(newp[3])) {
          JSINFO << MAGENTA << "THIS IS THE SECOND DAUGHTER";
          JSINFO << MAGENTA << plong << " " << s_t << " " << c_p << " "
                 << k_perp1[1];
          JSINFO << MAGENTA << newp[0] << " " << newp[1] << " " << newp[2]
                 << " " << newp[3];
          cin >> blurb;
        }

        VERBOSE(8) << MAGENTA << " D1 px = " << newp[1] << " py = " << newp[2]
                   << " pz = " << newp[3] << " E = " << newp[0];

        //newx[0] = time + deltaT;
        newx[0] = time;
        for (int j = 1; j <= 3; j++) {
          //newx[j] = pIn[i].x_in().comp(j) + (time + deltaT - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
          newx[j] = pIn[i].x_in().comp(j) +
                    (time - pIn[i].x_in().comp(0)) * velocity[j];
        }

        if (iSplit != 3) // not a photon
        {

          VERBOSE(8) << BOLDRED << " PiD - b = " << pid_b;
          pOut.push_back(Parton(0, pid_b, jet_stat, newp, newx));
          iout = pOut.size() - 1;
          pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
          pOut[iout].set_mean_form_time();
          ft = generate_L(pOut[iout].mean_form_time());
          pOut[iout].set_form_time(ft);
          pOut[iout].set_color(d2_col);
          pOut[iout].set_anti_color(d2_acol);
          pOut[iout].set_max_color(max_color);
          pOut[iout].set_min_color(pIn[i].min_color());
          pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

        } else // is a photon
        {
          VERBOSE(8) << BOLDRED << " is a photon PiD - b = " << pid_b;
          pOut.push_back(Photon(0, pid_b, jet_stat, newp, newx));
          iout = pOut.size() - 1;
          pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
          pOut[iout].set_mean_form_time();
          ft = generate_L(pOut[iout].mean_form_time());
          pOut[iout].set_form_time(1.0 / rounding_error);
          pOut[iout].set_color(0);
          pOut[iout].set_anti_color(0);
          pOut[iout].set_max_color(max_color);
          pOut[iout].set_min_color(pIn[i].min_color());
          pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

          //JSINFO << BOLDYELLOW << " A photon was made in Matter with px = " << pOut[iout].px() << " and sent to the framework " ;
        }

        VERBOSE(8) << BOLDRED << " virtuality of D 2 = " << pOut[iout].t();

      } else { // not time to split yet broadening it



        if (broadening_on) {

          double now_zeta =
              ((time + initRdotV + (time - initR0)) / std::sqrt(2)) *
              fmToGeVinv;
          qhat = fncQhat(now_zeta);
          if (now_temp > 0.1) {
            ehat = 0.0 * qhat / 4.0 / now_temp;
          } else {
            ehat = 0.0 * qhat / 4.0 / 0.3;
          }

          VERBOSE(8) << BOLDRED << " between splits broadening qhat = " << qhat
                     << " ehat = " << ehat << " and delT = " << delT;
          VERBOSE(8) << BOLDBLUE << " zeta at formation = " << zeta
                     << " zeta now = " << now_zeta;

          if ((!recoil_on) && (qhat > 0.0)) {
            double kt = 0;
            if (pIn[i].pid() == 21) // particle is a gluon
            {
              kt = generate_kt(qhat * 1.414 / 0.197, delT);
            } else if ((pIn[i].pid() < 6) &&
                       (pIn[i].pid() > -6)) // particle is a quark
            {
              kt = generate_kt(qhat * 1.414 / 0.197 * Cf / Ca,
                               delT); // scale down q-hat by Cf/Ca
            } else // a photon, or something else that does not have color
            {
              kt = 0;
            }

            JSDEBUG << " kt generated = " << kt
                    << " for qhat = " << qhat * 1.414 / 0.197
                    << " and delT = " << delT;

            double ktx, kty, ktz;
            ktx = kty = ktz = 0.0;
            double vx = initVx;
            double vy = initVy;
            double vz = initVz;

            bool solved = false;
            int trials = 0;
            while ((!solved) && (trials < 1000)) {

              if ((abs(vy) > approx) || (abs(vz) > approx)) {
                ktx =
                    kt * (1 - 2 * ZeroOneDistribution(*GetMt19937Generator()));

                JSDEBUG << " vx = " << vx << " vy = " << vy << " vz = " << vz;

                double rad = sqrt(
                    4 * ktx * ktx * vx * vx * vy * vy -
                    4 * (vy * vy + vz * vz) *
                        (ktx * ktx * (vx * vx + vz * vz) - kt * kt * vz * vz));

                double sol1 =
                    (-2 * ktx * vx * vy + rad) / 2 / (vy * vy + vz * vz);
                double sol2 =
                    (-2 * ktx * vx * vy - rad) / 2 / (vy * vy + vz * vz);

                kty = sol1;

                if ((ktx * ktx + sol1 * sol1) > kt * kt)
                  kty = sol2;

                if ((ktx * ktx + kty * kty) < kt * kt) {
                  ktz = sqrt(kt * kt - ktx * ktx - kty * kty);

                  double sign = ZeroOneDistribution(*GetMt19937Generator());

                  if (sign > 0.5)
                    ktz = -1 * ktz;

                  if (vz != 0)
                    ktz = (-1 * ktx * vx - kty * vy) / vz;

                  solved = true;
                }

              } else {
                ktx = 0;
                kty =
                    kt * (1 - 2 * ZeroOneDistribution(*GetMt19937Generator()));
                double sign = ZeroOneDistribution(*GetMt19937Generator());
                if (sign > 0.5)
                  kty = -1 * kty;
                ktz = sqrt(kt * kt - kty * kty);
                sign = ZeroOneDistribution(*GetMt19937Generator());
                if (sign > 0.5)
                  ktz = -1 * ktz;
              }
              trials++;
            }

            VERBOSE(8) << MAGENTA << " ktx = " << ktx << " kty = " << kty
                       << " ktz = " << ktz;

            double px = pIn[i].px();
            double py = pIn[i].py();
            double pz = pIn[i].pz();
            double energy = pIn[i].e();

            double p = sqrt(px * px + py * py + pz * pz);

            VERBOSE(8) << BOLDBLUE << " p before b & d, E = " << energy
                       << " pz = " << pz << " px = " << px << " py = " << py;

            px += ktx;
            py += kty;
            pz += ktz;

            double np = sqrt(px * px + py * py + pz * pz);
            JSDEBUG << " p = " << p << " np = " << np;

            double correction = 0.0;

            if (np * np > p * p)
              correction = sqrt(np * np - p * p);

            px -= correction * vx;
            py -= correction * vy;
            pz -= correction * vz;

            //                      double nnp = sqrt(px*px + py*py + pz*pz);

            /*                      if (abs(nnp-p)>abs(np-p))
                      {

                         VERBOSE(8) << MAGENTA << " negative condition invoked ! " ;

                          px += 2*(np-p)*vx;
                          py += 2*(np-p)*vy;
                          pz += 2*(np-p)*vz;
                      }
*/
            np = sqrt(px * px + py * py + pz * pz);
            JSDEBUG << "FINAL p = " << p << " np = " << np;

            pOut.push_back(pIn[i]);
            int iout = pOut.size() - 1;

            pOut[iout].reset_p(px, py, pz);

            double drag = ehat * delT;

            VERBOSE(8) << MAGENTA << " drag = " << drag
                       << " temperature = " << now_temp;

            if ((np > drag) && (energy > drag) &&
                (energy > sqrt(px * px + py * py + pz * pz))) {
              px -= drag * vx;
              py -= drag * vy;
              pz -= drag * vz;
              energy -= drag;
              pOut[iout].reset_momentum(px, py, pz, energy);
            }

            VERBOSE(8) << BOLDYELLOW << " p after b & d, E = " << energy
                       << " pz = " << pz << " px = " << px << " py = " << py;
          }

        } // end if(broadening_on)
          pOut.push_back(pIn[i]);
      }
    } 

    else { // virtuality too low lets broaden it



      if (broadening_on) {

        double now_zeta =
            ((time + initRdotV + (time - initR0)) / std::sqrt(2)) * fmToGeVinv;
        //double now_zeta = ( ( time + initRdotV )/std::sqrt(2) )*fmToGeVinv;
        qhat = fncQhat(now_zeta);
        if (now_temp > 0.1) {
          ehat = 0.0 * qhat / 4.0 / now_temp;
        } else {
          ehat = 0.0 * qhat / 4.0 / 0.3;
        }

        VERBOSE(8) << BOLDRED << " after splits broadening qhat = " << qhat
                   << " ehat = " << ehat << " and delT = " << delT;
        VERBOSE(8) << BOLDBLUE << " zeta at formation = " << zeta
                   << " zeta now = " << now_zeta;

        //JSINFO << " broadening qhat = " << qhat << " and delT = " << delT ;

        if ((!recoil_on) && (qhat > 0.0)) {
          double kt = 0.0;

          if (pIn[i].pid() == 21) {
            kt = generate_kt(qhat * 1.414 / 0.197, delT);
          } else {
            kt = generate_kt(qhat * 1.414 / 0.197 * Cf / Ca, delT);
          }

          JSDEBUG << " kt generated = " << kt
                  << " for qhat = " << qhat * 1.414 / 0.197
                  << " and delT = " << delT;

          double ktx, kty, ktz;
          ktx = kty = ktz = 0;
          double vx = initVx;
          double vy = initVy;
          double vz = initVz;

          bool solved = false;

          int trials = 0;
          while ((!solved) && (trials < 1000)) {
            if ((abs(vy) > approx) || (abs(vz) > approx)) {

              ktx = kt * (1 - 2 * ZeroOneDistribution(*GetMt19937Generator()));

              JSDEBUG << " vx = " << vx << " vy = " << vy << " vz = " << vz;

              double rad = sqrt(
                  4 * ktx * ktx * vx * vx * vy * vy -
                  4 * (vy * vy + vz * vz) *
                      (ktx * ktx * (vx * vx + vz * vz) - kt * kt * vz * vz));

              double sol1 =
                  (-2 * ktx * vx * vy + rad) / 2 / (vy * vy + vz * vz);
              double sol2 =
                  (-2 * ktx * vx * vy - rad) / 2 / (vy * vy + vz * vz);

              kty = sol1;

              if ((ktx * ktx + sol1 * sol1) > kt * kt)
                kty = sol2;

              if ((ktx * ktx + kty * kty) < kt * kt) {
                ktz = sqrt(kt * kt - ktx * ktx - kty * kty);

                double sign = ZeroOneDistribution(*GetMt19937Generator());
                if (sign > 0.5)
                  ktz = -1 * ktz;

                if (vz != 0)
                  ktz = (-1 * ktx * vx - kty * vy) / vz;

                solved = true;
              }
            } else {
              ktx = 0;
              kty = kt * (1 - 2 * ZeroOneDistribution(*GetMt19937Generator()));
              double sign = ZeroOneDistribution(*GetMt19937Generator());
              if (sign > 0.5)
                kty = -1 * kty;
              ktz = sqrt(kt * kt - kty * kty);
              sign = ZeroOneDistribution(*GetMt19937Generator());
              if (sign > 0.5)
                ktz = -1 * ktz;
            }
            trials++;
          }

          VERBOSE(8) << MAGENTA << " ktx = " << ktx << " kty = " << kty
                     << " ktz = " << ktz;

          double px = pIn[i].px();
          double py = pIn[i].py();
          double pz = pIn[i].pz();
          double energy = pIn[i].e();

          double p = sqrt(px * px + py * py + pz * pz);

          VERBOSE(8) << BOLDBLUE << " p before b & d, E = " << energy
                     << " pz = " << pz << " px = " << px << " py = " << py;

          px += ktx;
          py += kty;
          pz += ktz;

          double np = sqrt(px * px + py * py + pz * pz);
          JSDEBUG << " p = " << p << " np = " << np;
          double correction = 0.0;
          if (np * np > p * p)
            correction = sqrt(np * np - p * p);

          px -= correction * vx;
          py -= correction * vy;
          pz -= correction * vz;

          //  double nnp = sqrt(px*px + py*py + pz*pz);

          /*                  if (abs(nnp-p)>abs(np-p))
                  {
                      px += 2*(np-p)*vx;
                      py += 2*(np-p)*vy;
                      pz += 2*(np-p)*vz;
                  }
*/
          np = sqrt(px * px + py * py + pz * pz);
          JSDEBUG << "FINAL p = " << p << " nnnp = " << np;

          pOut.push_back(pIn[i]);
          int iout = pOut.size() - 1;

          pOut[iout].reset_p(px, py, pz);
          np = sqrt(px * px + py * py + pz * pz);
          double drag = ehat * delT;

          VERBOSE(8) << MAGENTA << " drag = " << drag
                     << " temperature = " << now_temp;

          if ((np > drag) && (energy > drag) &&
              (energy > sqrt(px * px + py * py + pz * pz))) {
            px -= drag * vx;
            py -= drag * vy;
            pz -= drag * vz;
            energy -= drag;
            pOut[iout].reset_momentum(px, py, pz, energy);
          }

          VERBOSE(8) << BOLDYELLOW << " p after b & d, E = " << energy
                     << " pz = " << pz << " px = " << px << " py = " << py;
        }

        //pOut.push_back(pIn[i]);
      }
    }

  } // particle loop
}

double Matter::generate_kt(double local_qhat, double dzeta) {
  double width = local_qhat * dzeta;
  //JSINFO<<"generate_kt local_qhat "<<local_qhat<<" dzeta "<<dzeta;
  double x, r;

  r = ZeroOneDistribution(*GetMt19937Generator());
  x = -log(1.0 - r);
  if (x < 0.0)
    throw std::runtime_error(" k_t^2 < 0.0 ");
  double kt = sqrt(x * local_qhat * dzeta);

  return (kt);
}

double Matter::generate_angle() {
  double ang, r, blurb;

  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //    r = mtrand1();

  //    cout << " r = " << r << endl;
  //    cin >> blurb ;

  ang = r * 2.0 * pi;

  return (ang);
}

double Matter::generate_vac_t(int p_id, double nu, double t0, double t,
                              double loc_a, int is) {
  double r, z, ratio, diff, scale, t_low, t_hi, t_mid, numer, denom, test;
  //JSINFO<<"generate_vac_t pid "<<p_id<<" nu "<<nu<<" t0 "<<t0<<" t "<<t<<" loc_a "<<loc_a<<" is "<<is;
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //        r = mtrand1();

  if ((r >= 1.0) || (r <= 0.0)) {
    throw std::runtime_error(
        "error in random number in t *GetMt19937Generator()");
  }

  ratio = 1.0;

  diff = (ratio - r) / r;

  t_low = 2.0 * t0;

  t_hi = t;

  tscale = t; //for virtuality dependent q-hat
  //    cout << " in gen_vac_t : t_low , t_hi = " << t_low << "  " << t_hi << endl;
  //    cin >> test ;

  if (p_id == gid) {
    numer = SudakovFunction::sudakov_Pgg(t0, t, loc_a, nu,*this) *
            std::pow(SudakovFunction::sudakov_Pqq(t0, t, loc_a, nu,*this), nf);

    if ((is != 1) &&
        (is !=
         2)) // there is almost no use of is = 2, the `is' in this function is redundant, just for consistency
    {
      throw std::runtime_error(" error in isp ");
    }
  } else {
    if ((is != 0) &&
        (is !=
         3)) // there is almost no use of is = 3, the `is' in this function is redundant, just for consistency
    {
      throw std::runtime_error("error in isp in quark split");
    }
    double relative_charge = 0.0;
    if ((std::abs(p_id) > 0) && (std::abs(p_id) < 4))
      relative_charge = 1.0 / 9.0;
    if (std::abs(p_id) == 2)
      relative_charge = relative_charge * 4.0;

    numer = SudakovFunction::sudakov_Pqg(t0, t, loc_a, nu,*this) *
            std::pow(SudakovFunction::sudakov_Pqp(t0, t, loc_a, nu,*this), relative_charge);
  }

  t_mid = t_low;

  if (numer > r) {
    // cout << " numer > r, i.e. ; " << numer << " > " << r << endl ;
    return (t_mid);
  }

  //t_mid = (t_low+t_hi)/2.0 ;

  scale = t0;

  //   cout << " s_approx, s_error = " << s_approx << "  " << s_error << endl;

  do {

    t_mid = (t_low + t_hi) / 2.0;

    if (p_id == gid) {
      denom = SudakovFunction::sudakov_Pgg(t0, t_mid, loc_a, nu,*this) *
              std::pow(SudakovFunction::sudakov_Pqq(t0, t_mid, loc_a, nu,*this), nf);
      if ((is != 1) && (is != 2)) {
        throw std::runtime_error(" error in isp numerator");
      }

    } else {
      if ((is != 0) && (is != 3)) {
        throw std::runtime_error(" error in isp in quark split numerator  ");
      }
      double relative_charge = 0.0;
      if ((std::abs(p_id) > 0) && (std::abs(p_id) < 4))
        relative_charge = 1.0 / 9.0;
      if (std::abs(p_id) == 2)
        relative_charge = relative_charge * 4.0;

      denom = SudakovFunction::sudakov_Pqg(t0, t_mid, loc_a, nu,*this) *
              std::pow(SudakovFunction::sudakov_Pqp(t0, t_mid, loc_a, nu,*this), relative_charge);
    }

    ratio = numer / denom;

    diff = (ratio - r) / r;

    //       cout << "num, den, r = " << numer << " "<< denom << " " << r << " " << endl;
    //       cout << "diff, t_mid = " << diff << " " << t_mid << endl;
    //       cout << " t_low, t_hi = " << t_low << "  " << t_hi << endl;
    //       cin >> test ;

    if (diff < 0.0) {
      t_low = t_mid;
      //t_mid = (t_low + t_hi)/2.0;
    } else {
      t_hi = t_mid;
      //t_mid = (t_low + t_hi)/2.0;
    }

  } while ((abs(diff) > s_approx) && (abs(t_hi - t_low) / t_hi > s_error));

  return (t_mid);
}

double Matter::generate_vac_t_w_M(int p_id, double M, double nu, double t0,
                                  double t, double loc_a, int is) {
  double r, z, ratio, diff, scale, t_low_M0, t_low_MM, t_low_00, t_hi_M0,
      t_hi_MM, t_hi_00, t_mid_M0, t_mid_MM, t_mid_00, numer, denom, test;
  //JSINFO<<"generate_vac_t_w_M pid "<<p_id<<" nu "<<nu<<" t0 "<<t0<<" t "<<t<<" loc_a "<<loc_a<<" is "<<is;
  double M_charm = PythiaFunction.particleData.m0(cid);
  double M_bottom = PythiaFunction.particleData.m0(bid);

  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //        r = mtrand1();

  if ((r >= 1.0) || (r <= 0.0)) {
    throw std::runtime_error(
        "error in random number in t *GetMt19937Generator()");
  }

  ratio = 1.0;

  diff = (ratio - r) / r;

  //    if (t0<M*M) t0 = M*M;

  t_low_M0 = t0 * (1.0 + std::sqrt(1.0 + 2.0 * M * M / t0));
  t_low_MM = 2.0 * (M * M + t0);
  t_low_00 = 2.0 * t0;

  t_hi_M0 = t;
  t_hi_MM = t;
  t_hi_00 = t;

  tscale = t; //for virtuality dependent q-hat

  VERBOSE(1) << MAGENTA << " in gen_vac_t_w_M : t_low , t_hi = " << t_low_M0 << "  " << t ;
  //cin >> test ;

  if (p_id == gid) {
    if (t < 2.0 * (M_charm * M_charm + t0)) {
      numer = SudakovFunction::sudakov_Pgg(t0, t, loc_a, nu,*this) *
              std::pow(SudakovFunction::sudakov_Pqq(t0, t, loc_a, nu,*this), 3.0);
    } else if (t < 2.0 * (M_bottom * M_bottom + t0)) {
      numer = SudakovFunction::sudakov_Pgg(t0, t, loc_a, nu,*this) *
              std::pow(SudakovFunction::sudakov_Pqq(t0, t, loc_a, nu,*this), 3.0) *
              SudakovFunction::sudakov_Pqq_w_M_vac_only(M_charm, t0, t, loc_a, nu,*this);
    } else {
      numer = SudakovFunction::sudakov_Pgg(t0, t, loc_a, nu,*this) *
              std::pow(SudakovFunction::sudakov_Pqq(t0, t, loc_a, nu,*this), 3.0) *
              SudakovFunction::sudakov_Pqq_w_M_vac_only(M_charm, t0, t, loc_a, nu,*this) *
              SudakovFunction::sudakov_Pqq_w_M_vac_only(M_bottom, t0, t, loc_a, nu,*this);
    }
    // JSINFO << BOLDYELLOW << " numer = " << numer ;

    // double blurb;
    // cin >> blurb;

    if ((is != 1) && (is != 2) && (is != 4) && (is != 5)) {
      throw std::runtime_error(" error in isp ");
    }
  } else {
    if (is != 0) {
      throw std::runtime_error("error in isp in quark split");
    }
    if (((int)std::abs((double)p_id)) == 4 ||
        ((int)std::abs((double)p_id)) == 5) {
      numer = SudakovFunction::sudakov_Pqg_w_M(M, t0, t, loc_a, nu,*this);

      //    std::ofstream Sud_dist;
      //    Sud_dist.open("Sud_dist.dat", std::ios::app);
      //    Sud_dist << t << "  " << numer << "  " << sudakov_Pqg(t0,t,loc_a,nu) << "  " << r << " t_low_MO = " << t_low_M0 << endl;
      //    Sud_dist.close();
    } else {
      numer = SudakovFunction::sudakov_Pqg(t0, t, loc_a, nu,*this);
    }
  }

  t_mid_M0 = t_low_M0;
  t_mid_MM = t_low_MM;
  t_mid_00 = t_low_00;

  if (numer > r) {
    // cout << " numer > r, i.e. ; " << numer << " > " << r << endl ;
    if (std::fabs(p_id) == cid || std::fabs(p_id) == bid)
      return (t_mid_M0);
    else
      return (t_mid_00);
  }

  //t_mid = (t_low+t_hi)/2.0 ;

  scale = t0;

  //   cout << " s_approx, s_error = " << s_approx << "  " << s_error << endl;

  bool exit_condition = true;

  do {
    t_mid_M0 = (t_low_M0 + t_hi_M0) / 2.0;
    t_mid_MM = (t_low_MM + t_hi_MM) / 2.0;
    t_mid_00 = (t_low_00 + t_hi_00) / 2.0;

    if (p_id == gid) {
      if (t_mid_00 < 2.0 * (M_charm * M_charm + t0)) {
        denom = SudakovFunction::sudakov_Pgg(t0, t_mid_00, loc_a, nu,*this) *
                std::pow(SudakovFunction::sudakov_Pqq(t0, t_mid_00, loc_a, nu,*this), 3.0);
      } else if (t_mid_00 < 2.0 * (M_bottom * M_bottom + t0)) {
        denom = SudakovFunction::sudakov_Pgg(t0, t_mid_00, loc_a, nu,*this) *
                std::pow(SudakovFunction::sudakov_Pqq(t0, t_mid_00, loc_a, nu,*this), 3.0) *
                SudakovFunction::sudakov_Pqq_w_M_vac_only(M_charm, t0, t_mid_00, loc_a, nu,*this);
      } else {
        denom = SudakovFunction::sudakov_Pgg(t0, t_mid_00, loc_a, nu,*this) *
                std::pow(SudakovFunction::sudakov_Pqq(t0, t_mid_00, loc_a, nu,*this), 3.0) *
                SudakovFunction::sudakov_Pqq_w_M_vac_only(M_charm, t0, t_mid_00, loc_a, nu,*this) *
                SudakovFunction::sudakov_Pqq_w_M_vac_only(M_bottom, t0, t_mid_00, loc_a, nu,*this);
      }
      if ((is != 1) && (is != 2) && (is != 4) && (is != 5)) {
        throw std::runtime_error(" error in isp numerator");
      }
    } else {
      if (is != 0) {
        throw std::runtime_error(" error in isp in quark split numerator  ");
      }
      if (((int)std::abs((double)p_id)) == 4 ||
          ((int)std::abs((double)p_id)) == 5) {
        VERBOSE(8) << BOLDYELLOW << " In generate_vac_t_w_M,  t0 = " << t0
                   << " t_mid_M0 = " << t_mid_M0;
        denom = SudakovFunction::sudakov_Pqg_w_M(M, t0, t_mid_M0, loc_a, nu,*this);
      } else {
        denom = SudakovFunction::sudakov_Pqg(t0, t_mid_00, loc_a, nu,*this);
      }
    }

    ratio = numer / denom;

    diff = (ratio - r) / r;

    if (diff < 0.0) {
      t_low_M0 = t_mid_M0;
      t_low_MM = t_mid_MM;
      t_low_00 = t_mid_00;

      //t_mid = (t_low + t_hi)/2.0;
    } else {
      t_hi_M0 = t_mid_M0;
      t_hi_MM = t_mid_MM;
      t_hi_00 = t_mid_00;

      //t_mid = (t_low + t_hi)/2.0;
    }

    if (std::abs(p_id) == cid || std::abs(p_id) == bid)
      exit_condition = (std::abs(diff) < s_approx) ||
                       (std::abs(t_hi_M0 - t_low_M0) / t_hi_M0 < s_error);
    if (p_id == gid)
      exit_condition = (std::abs(diff) < s_approx) ||
                       (std::abs(t_hi_00 - t_low_00) / t_hi_00 < s_error);
    // need to think about the second statement in the gluon exit condition.

  } while (!exit_condition);

  if (std::fabs(p_id) == cid || std::fabs(p_id) == bid)
    return (t_mid_M0);
  else
    return (t_mid_00);
}

/*



  New function



*/

double Matter::generate_vac_z(int p_id, double t0, double t, double loc_b,
                              double nu, int is) {
  double r, z, ratio, diff, e, numer1, numer2, numer, denom, z_low, z_hi, z_mid,
      test;
  //JSINFO<<"generate_vac_z pid "<<p_id<<" nu "<<nu<<" t0 "<<t0<<" t "<<t<<" loc_b "<<loc_b<<" is "<<is;

  r = ZeroOneDistribution(*GetMt19937Generator());

  if ((r > 1) || (r < 0)) {
    throw std::runtime_error(
        " error in random number in z *GetMt19937Generator()");
  }

  ratio = 1.0;

  diff = (ratio - r) / r;

  e = t0 / t;

  if (e > 0.5) {
    throw std::runtime_error(" error in epsilon");
  }

  z_low = e;

  z_hi = double(1.0) - e;

  if (p_id == gid) {
    if (is == 1) {
      denom = SplittingFunction::P_z_gg_int(z_low, z_hi, loc_b, t, 2.0 * nu / t, nu,*this);
      //JSINFO<<" "<<denom;
    } else {
      denom = SplittingFunction::P_z_qq_int(z_low, z_hi, loc_b, t, 2.0 * nu / t, nu,*this);
    }

  } else if ((p_id != gid) && (is == 0)) {
    denom = SplittingFunction::P_z_qg_int(z_low, z_hi, loc_b, t, 2.0 * nu / t, nu,*this);
  } else if ((p_id != gid) && (is == 3)) {
    denom = SplittingFunction::P_z_qp_int(z_low, z_hi, loc_b, t, 2.0 * nu / t, nu,*this);
  } else {
    throw std::runtime_error(
        " I do not understand your combination of particle id and split id in "
        "denominator of generate_z ");
  }

  //z_mid = (z_low + z_hi)/2.0 ;

  int itcounter = 0;
  // cout << " generate_vac_z called with p_id = " << p_id << " t0 = " << t0 << " t = " << t << " loc_b=" << loc_b<< " nu = " <<  nu << " is = " << is << endl;

  do { // Getting stuck in here for some reason
    if (itcounter++ > 10000) {
      cout << " in here"
           << " abs(diff) = " << abs(diff) << "  approx = " << approx
           << "  r = " << r << "  zmid = " << z_mid << "  denom = " << denom
           << "  numer = " << numer << "  e = " << e << "   " << numer / denom
           << endl;
      throw std::runtime_error("Stuck in endless loop");
    }

    z_mid = (z_low + z_hi) / 2.0;

    if (p_id == gid) {
      if (is == 1) {
        numer = SplittingFunction::P_z_gg_int(e, z_mid, loc_b, t, 2.0 * nu / t, nu,*this);
        //JSINFO<<" "<<numer;
      } else {
        numer = SplittingFunction::P_z_qq_int(e, z_mid, loc_b, t, 2.0 * nu / t, nu,*this);
      }
    } else if ((p_id != gid) && (is == 0)) {
      numer = SplittingFunction::P_z_qg_int(e, z_mid, loc_b, t, 2.0 * nu / t, nu,*this);
    } else if ((p_id != gid) && (is == 3)) {
      numer = SplittingFunction::P_z_qp_int(e, z_mid, loc_b, t, 2.0 * nu / t, nu,*this);
    } else {
      throw std::runtime_error(
          " I do not understand your combination of particle id and split id "
          "in numerator of generate_z ");
    }
    ratio = numer / denom;
    diff = (ratio - r) / r;

    if (diff > 0.0) {
      z_hi = z_mid;
      //z_mid = (z_low + z_hi)/2.0;
    } else {
      z_low = z_mid;
      //z_mid = (z_low + z_hi)/2.0 ;
    }

  } while ((abs(diff) > approx) && (abs(z_hi - z_low) / z_hi > s_error));

  return (z_mid);
}

double Matter::generate_vac_z_w_M(int p_id, double M, double t0, double t,
                                  double loc_b, double nu, int is) {
  double r, z, ratio, diff, e, numer1, numer2, numer, denom, z_low_00, z_low_M0,
      z_low_MM, z_hi_00, z_hi_M0, z_hi_MM, z_mid_00, z_mid_M0, z_mid_MM, test,
      z_min_00, z_min_M0, z_min_MM;
  //JSINFO<<"generate_vac_z_w_M pid "<<p_id<<" nu "<<nu<<" t0 "<<t0<<" t "<<t<<" loc_b "<<loc_b<<" is "<<is;
  r = ZeroOneDistribution(*GetMt19937Generator());

  //    if (t0<M*M) t0 = M*M;

  if ((r > 1) || (r < 0)) {
    throw std::runtime_error(
        " error in random number in z *GetMt19937Generator()");
  }

  ratio = 1.0;

  diff = (ratio - r) / r;

  e = t0 / t;

  double M2 = M * M;

  if (e > 0.5) {
    throw std::runtime_error(" error in epsilon");
  }

  z_low_M0 = e + M2 / (t + M2);
  z_low_00 = e;
  z_low_MM = e + M2 / t;

  z_hi_00 = double(1.0) - e;
  z_hi_M0 = double(1.0) - e;
  z_hi_MM = double(1.0) - e - M2 / t;

  z_min_00 = z_low_00;
  z_min_M0 = z_low_M0;
  z_min_MM = z_low_MM;

  // JSINFO << BOLDYELLOW << " r = " << r << " 1st z_low_00 = " << z_low_00 << "1st z_hi_00 = " << z_hi_00 << " M = " << M << " is = " << is << " pid = " << p_id;

  if (p_id == gid) {
    if (is == 1) // for a gluon to 2 gluons
    {
      denom = SplittingFunction::P_z_gg_int(z_min_00, z_hi_00, loc_b, t, 2.0 * nu / t, nu,*this);
      //JSINFO<<" "<<denom;
      //JSINFO << MAGENTA << " denom = " << denom ;
    } else {
      if (is > 3) //THIS IS FOR GLUON-> HEAVY Q-QBAR
      {
        denom = SplittingFunction::P_z_qq_int_w_M_vac_only(M, z_min_MM, z_hi_MM, loc_b, t,
                                        2.0 * nu / t, nu,*this);
      } else if (is == 2) // For a gluon to light quark anti-quark
      {
        denom = SplittingFunction::P_z_qq_int(z_min_00, z_hi_00, loc_b, t, 2.0 * nu / t, nu,*this);
      }
    }
  } else {
    if ((std::abs(p_id) == 4) || (std::abs(p_id) == 5)) {
      denom = SplittingFunction::P_z_qg_int_w_M(M, z_min_M0, z_hi_M0, loc_b, t, 2.0 * nu / t, nu,*this);
    } else {
      denom = SplittingFunction::P_z_qg_int(z_min_00, z_hi_00, loc_b, t, 2.0 * nu / t, nu,*this);
    }
  }

  //z_mid = (z_low + z_hi)/2.0 ;

  int itcounter = 0;
  //JSINFO << BOLDYELLOW << " generate_vac_z called with p_id = " << p_id << " t0 = " << t0 << " t = " << t << " loc_b=" << loc_b<< " nu = " <<  nu << " is = " << is ;

  do { // Getting stuck in here for some reason

    //JSINFO << BOLDYELLOW << " r = " << r << " z_low_00 in loop = " << z_low_00 << " z_hi_00 in loop = " << z_hi_00 << " M = " << M << " is = " << is << " pid = " << p_id;

    if (itcounter++ > 10000) {
      cout << " in here"
           << " abs(diff) = " << abs(diff) << "  approx = " << approx
           << "  r = " << r << "  zmid = " << z_mid_M0 << "  denom = " << denom
           << "  numer = " << numer << "  e = " << e << "   " << numer / denom
           << endl;
      throw std::runtime_error("Stuck in endless loop");
    }

    z_mid_00 = (z_low_00 + z_hi_00) / 2.0;
    z_mid_M0 = (z_low_M0 + z_hi_M0) / 2.0;
    z_mid_MM = (z_low_MM + z_hi_MM) / 2.0;

    if (p_id == gid) {
      if (is == 1) {
        numer = SplittingFunction::P_z_gg_int(z_min_00, z_mid_00, loc_b, t, 2.0 * nu / t, nu,*this);
        //JSINFO<<" "<<numer;
        //JSINFO << MAGENTA << " numer = " << numer ;
      } else {
        if (is >
            3) // 3 is quark radiating a photon, 4,5 is gluon radiating a c cbar, b bbar
        {
          numer = SplittingFunction::P_z_qq_int_w_M_vac_only(M, z_min_MM, z_mid_MM, loc_b, t,
                                          2.0 * nu / t, nu,*this);
        } else if (is == 2) // gluon to light quark anti-quark
        {
          numer = SplittingFunction::P_z_qq_int(z_min_00, z_mid_00, loc_b, t, 2.0 * nu / t, nu,*this);
        }
      }
    } else {

      if ((std::abs(p_id) == 4) || (std::abs(p_id) == 5)) {
        numer =
            SplittingFunction::P_z_qg_int_w_M(M, z_min_M0, z_mid_M0, loc_b, t, 2.0 * nu / t, nu,*this);
      } else {
        numer = SplittingFunction::P_z_qg_int(z_min_00, z_mid_00, loc_b, t, 2.0 * nu / t, nu, *this);
      }
    }
    ratio = numer / denom;
    diff = (ratio - r) / r;
    //JSINFO<< BOLDYELLOW << "num = " << numer << " denom = "<< denom << " ratio = " << numer/denom << " r = " << r ;
    // JSINFO << BOLDYELLOW << " diff = " << diff << " z_mid = " << z_mid_00 ;
    //	cin >> test ;

    if ((diff == 0.0) && ((p_id == gid) || (std::abs(p_id) < 4)))
      return (z_mid_00);
    if ((diff == 0.0) && ((std::abs(p_id) == 4) || (std::abs(p_id) == 5)))
      return (z_mid_M0);

    if (diff > 0.0) {
      z_hi_M0 = z_mid_M0;
      z_hi_00 = z_mid_00;
      z_hi_MM = z_mid_MM;
      //z_mid = (z_low + z_hi)/2.0;
    } else {
      z_low_M0 = z_mid_M0;
      z_low_00 = z_mid_00;
      z_low_MM = z_mid_MM;
      //z_mid = (z_low + z_hi)/2.0 ;
    }

  } while (
      ((std::abs(p_id) == cid || std::abs(p_id) == bid) &&
       (abs(diff) > approx) && (abs(z_hi_M0 - z_low_M0) / z_hi_M0 > s_error)) ||
      ((std::abs(p_id) == uid || std::abs(p_id) == did ||
        std::abs(p_id) == sid) &&
       (abs(diff) > approx) && (abs(z_hi_00 - z_low_00) / z_hi_00 > s_error)) ||
      ((std::abs(p_id) == gid && is >= 1 && is < 3) && (abs(diff) > approx) &&
       (abs(z_hi_00 - z_low_00) / z_hi_00 > s_error)) ||
      ((std::abs(p_id) == gid && is > 3) && (abs(diff) > approx) &&
       (abs(z_hi_MM - z_low_MM) / z_hi_MM > s_error)));

  // JSINFO << BOLDRED << " z_mid_00 = " << z_mid_00 ;

  if (p_id == gid && (is == 1 || is == 2))
    return (z_mid_00);
  else if (p_id == gid && (is == 4 || is == 5))
    return (z_mid_MM);
  else if (std::abs(p_id) == uid || std::abs(p_id) == did ||
           std::abs(p_id) == sid)
    return (z_mid_00);
  else
    return (z_mid_M0);
}

double Matter::generate_L(double form_time) {
  double r, x_low, x_high, x, diff, span, val, arg, norm;
  //JSINFO<<"generate_L form_time "<<form_time;
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //    r = mtrand1();

  if ((r > 1) || (r < 0)) {
    throw std::runtime_error(
        " error in random number in z *GetMt19937Generator()");
  }

  x_low = 0;

  x_high = 8.0 * form_time;
  // the value of x_high is slightly arbitrary, the erf function is more or less zero at this distance.
  // picking 10*form_time will not lead to any different results

  x = (x_low + x_high) / 2.0;

  span = (x_high - x_low) / x_high;

  arg = x / form_time / std::sqrt(pi);

  val = std::erf(arg);

  diff = std::abs(val - r);

  while ((diff > approx) && (span > error)) {
    if ((val - r) > 0.0) {
      x_high = x;
    } else {
      x_low = x;
    }

    x = (x_low + x_high) / 2.0;

    arg = x / form_time / std::sqrt(pi);

    val = std::erf(arg);

    diff = std::abs(val - r);

    span = (x_high - x_low) / x_high;
  }

  //	cout << " random number for dist = " << r << " distance generated = " << x << endl;

  return (x);
}



double Matter::alpha_s(double q2) {
  double a, L2, q24, c_nf;

  L2 = std::pow(Lambda_QCD, 2);

  q24 = q2 / 4.0;

  c_nf = nf;

  if (q24 > 4.0) {
    c_nf = 4;
  }

  if (q24 > 64.0) {
    c_nf = 5;
  }

  if (q24 > L2) {
    a = 12.0 * pi / (11.0 * Nc - 2.0 * c_nf) / std::log(q24 / L2);
  } else {
    JSWARN << " alpha too large ";
    a = 0.6;
  }

  return (a);
}

double Matter::profile(double zeta) {
  double prof;

  /*  Modify the next set of lines to get a particular profile in brick test or further modify the profile for hydro to introduce coherence effects */

  prof = 1.0;

  return (prof);
}

////////////////////////////////////////////////////////////////////////////////////////

double Matter::fillQhatTab(double y) {

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
double Matter::GeneralQhatFunction(int QhatParametrization, double Temperature, double EntropyDensity, double FixAlphas, double Qhat0, double E, double muSquare)
{
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
double Matter::RunningAlphaS(double muSquare)
{
  int ActiveFlavor=3;
  double Square_Lambda_QCD_HTL = exp( -12.0*pi/( (33 - 2*ActiveFlavor)*alphas) );
  double ans = 12.0*pi/( (33.0- 2.0*ActiveFlavor)*log(muSquare/Square_Lambda_QCD_HTL) );
  if(muSquare < 1.0) {ans=alphas; }

  VERBOSE(8)<<"Fixed-alphaS="<<alphas<<", Lambda_QCD_HTL="<<sqrt(Square_Lambda_QCD_HTL)<<", mu2="<<muSquare<<", Running alpha_s"<<ans;
  return ans;
}
//////////////////////////////////////////////////////////////////////////////////////
///////////// Virtuality dependent prefactor for q-hat function /////////////////////////
double Matter::VirtualityQhatFunction(int QhatParametrization,  double enerLoc, double muSquare)
{
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
////////////Modification of elastic scattering probability due to modified q-hat//////
double Matter::ModifiedProbability(int QhatParametrization, double tempLoc, double sdLoc, double enerLoc, double muSquare)
{
  double ModifiedAlphas=0;double qhatLoc=0;
  double ScaleNet=2*enerLoc*tempLoc;
  if(ScaleNet <1.0) { ScaleNet=1.0; }

  switch(QhatParametrization)
    {
      //For HTL q-hat formula with all alpha_s as constant and controlled by XML
    case 0:
      ModifiedAlphas = alphas;
      break;

      //For HTL q-hat where one alpha_s at scale muS=2ET and second alpha_s at muS=DebyeMassSquare is fit parameter
    case 1:
      ModifiedAlphas = RunningAlphaS(ScaleNet);
      break;

      //Constant q-hat case, alphas is computed using HTL formula with all alpha_s fixed
    case 2:
      qhatLoc = qhat0*0.1973;
      ModifiedAlphas = solve_alphas(qhatLoc, enerLoc, tempLoc);
      break;

      //For q-hat goes as T^3, alphas is computed using HTL formula with all alpha_s fixed
    case 3:
      qhatLoc = qhat0*pow(tempLoc/0.3,3)*0.1973; // w.r.t T=0.3 GeV
      ModifiedAlphas = solve_alphas(qhatLoc, enerLoc, tempLoc);
      break;

      //For q-hat goes as entropy density
    case 4:
      qhatLoc = qhat0*(sdLoc/96.0)*0.1973; // w.r.t S0=96 fm^-3
      ModifiedAlphas = solve_alphas(qhatLoc, enerLoc, tempLoc);
      break;

      //For HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat
      //Function is 1 / (1+A*pow(log(Q^2),2)+B*pow(log(Q^2),4))
    case 5:
      ModifiedAlphas = RunningAlphaS(ScaleNet)*VirtualityQhatFunction(5,  enerLoc, muSquare) ;
      break;
      //HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat
      //Function is int^{1}_{xB} e^{-ax} / (1+A*pow(log(Q^2),1)+B*pow(log(Q^2),2))
    case 6:
      ModifiedAlphas = RunningAlphaS(ScaleNet)*VirtualityQhatFunction(6,  enerLoc, muSquare) ;
      break;

      //HTL q-hat multiplied by Virtuality dependent function to mimic PDF-Scale dependent q-hat
      //Function is int^{1}_{xB} x^{a}(1-x)^{b} / (1+A*pow(log(Q^2),1)+B*pow(log(Q^2),2))
    case 7:
      ModifiedAlphas = RunningAlphaS(ScaleNet)*VirtualityQhatFunction(7,  enerLoc, muSquare) ;
      break;

    default:
      JSINFO<<"q-hat Parametrization "<<QhatParametrization<<" is not used,Elastic scattering alphas will be set to zero";
    }
  //JSINFO<<"q-hat Parametrization "<<QhatParametrization<<" modified alphas="<<ModifiedAlphas;
  return ModifiedAlphas;
}

//x integration  of QGP-PDF from xB to 1
double Matter::IntegralPDF(double xB, double a, double b)
{
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


double Matter::fncQhat(double zeta) {
  if (in_vac)
    return (0.0);

  double tStep = 0.1;
  //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
  int indexZeta =
      (int)((sqrt(2.0) * zeta / 5.0 - initRdotV + initR0) / 2.0 / tStep +
            0.5); // zeta was in 1/GeV and light cone coordinate

  //if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;
  if (indexZeta >= dimQhatTab)
    return (0);

  double avrQhat = qhatTab1D[indexZeta]*VirtualityQhatFunction(QhatParametrizationType, initEner, tscale);
  return (avrQhat);
}

//////////////////////////////////////////////////////////////////////////////////////

double Matter::fncAvrQhat(double zeta, double tau) {

  if (in_vac)
    return (0.0);

  double tStep = 0.1;
  //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
  int indexZeta =
      (int)((sqrt(2.0) * zeta / 5.0 - initRdotV + initR0) / 2.0 / tStep +
            0.5); // zeta was in 1/GeV and light cone coordinate
  int indexTau = (int)(tau / sqrt(2.0) / 5.0 / tStep +
                       0.5); // tau was in 1/GeV and light cone coordinate

  // if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;
  if (indexZeta >= dimQhatTab)
    return (0);
  if (indexTau >= dimQhatTab)
    indexTau = dimQhatTab - 1;

  double avrQhat = qhatTab2D[indexZeta][indexTau]*VirtualityQhatFunction(QhatParametrizationType, initEner, tscale);
  return (avrQhat);
}

//////////////////////////////////////////////////////////////////////////////////////

void Matter::flavor(int &CT, int &KATT0, int &KATT2, int &KATT3,
                    unsigned int &max_color, unsigned int &color0,
                    unsigned int &anti_color0, unsigned int &color2,
                    unsigned int &anti_color2, unsigned int &color3,
                    unsigned int &anti_color3) {

  int vb[7] = {0};
  int b = 0;
  int KATT00 = KATT0;
  unsigned int backup_color0 = color0;
  unsigned int backup_anti_color0 = anti_color0;

  vb[1] = 1;
  vb[2] = 2;
  vb[3] = 3;
  vb[4] = -1;
  vb[5] = -2;
  vb[6] = -3;

  if (KATT00 == 21) { //.....for gluon
    double R1 = 16.0; // gg->gg DOF_g
    double R2 = 0.0;  // gg->qqbar don't consider this channel in Matter
    double R3 = 6.0 * 6 * 4 / 9; // gq->gq or gqbar->gqbar flavor*DOF_q*factor
    double R0 = R1 + R3;

    double a = ran0(&NUM1);

    if (a <= R1 / R0) { // gg->gg
      CT = 1;
      KATT3 = 21;
      KATT2 = 21;
      //KATT0=KATT0;
      max_color++;
      color0 = backup_color0;
      anti_color0 = max_color;
      max_color++;
      color2 = anti_color0;
      anti_color2 = max_color;
      color3 = backup_anti_color0;
      anti_color3 = max_color;
    } else { // gq->gq
      CT = 3;
      b = floor(ran0(&NUM1) * 6 + 1);
      if (b == 7)
        b = 6;
      KATT3 = vb[b];
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT3 > 0) { // gq->gq
        max_color++;
        color0 = backup_color0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { // gqbar->gqbar
        max_color++;
        color0 = max_color;
        anti_color0 = backup_anti_color0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      }
    }
  } else if (abs(KATT00) == 4) { // for charm quarks
    double R1 = 6.0 * 6 * 4 / 9; // Qq->Qq
    double R2 = 16.0;            // Qg->Qg DOF_ag
    double R00 = R1 + R2;

    double a = ran0(&NUM1);

    if (a <= R2 / R00) { // Qg->Qg
      CT = 12;
      KATT3 = 21;
      KATT2 = 21;
      if (KATT00 > 0) { // Qg->Qg
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        max_color++;
        color2 = max_color;
        anti_color2 = color0;
        color3 = max_color;
        anti_color3 = backup_color0;
      } else { // Qbarg->Qbarg
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        max_color++;
        color2 = anti_color0;
        anti_color2 = max_color;
        color3 = backup_anti_color0;
        anti_color3 = max_color;
      }
    } else { // Qq->Qq
      CT = 11;
      b = floor(ran0(&NUM1) * 6 + 1);
      if (b == 7)
        b = 6;
      KATT3 = vb[b];
      KATT2 = KATT3;
      if (KATT00 > 0 && KATT2 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else if (KATT00 > 0 && KATT2 < 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else if (KATT00 < 0 && KATT2 > 0) { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    }

  } else {                       //.....for quark and antiquark (light)
    double R3 = 16.0;            // qg->qg DOF_g
    double R4 = 4.0 * 6 * 4 / 9; // qq'->qq' scatter with other species
    double R5 = 1.0 * 6 * 4 / 9; // qq->qq scatter with itself
    double R6 =
        0.0; // qqbar->q'qbar' to other final state species, don't consider in Matter
    double R7 = 1.0 * 6 * 4 / 9; // qqbar->qqbar scatter with its anti-particle
    double R8 = 0.0;             // qqbar->gg don't consider in Matter
    double R00 = R3 + R4 + R5 + R7;

    double a = ran0(&NUM1);
    if (a <= R3 / R00) { // qg->qg
      CT = 13;
      KATT3 = 21;
      KATT2 = 21;
      //KATT0=KATT0;
      if (KATT00 > 0) { // qg->qg
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        max_color++;
        color2 = max_color;
        anti_color2 = color0;
        color3 = max_color;
        anti_color3 = backup_color0;
      } else { // qbarg->qbarg
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        max_color++;
        color2 = anti_color0;
        anti_color2 = max_color;
        color3 = backup_anti_color0;
        anti_color3 = max_color;
      }
    } else if (a <= (R3 + R4) / R00) { // qq'->qq'
      CT = 4;
      do {
        b = floor(ran0(&NUM1) * 6 + 1);
        if (b == 7)
          b = 6;
        KATT3 = vb[b];
      } while (KATT3 == KATT0 || KATT3 == -KATT0);
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT00 > 0 && KATT2 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else if (KATT00 > 0 && KATT2 < 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else if (KATT00 < 0 && KATT2 > 0) { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    } else if (a <= (R3 + R4 + R5) / R00) { // scatter with itself
      CT = 5;
      KATT3 = KATT0;
      KATT2 = KATT0;
      //KATT0=KATT0;
      if (KATT00 > 0) { // qq->qq
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = backup_color0;
        anti_color2 = 0;
        color3 = max_color;
        anti_color3 = 0;
      } else { //qbarqbar->qbarqbar
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = 0;
        anti_color2 = backup_anti_color0;
        color3 = 0;
        anti_color3 = max_color;
      }
    } else { // scatter with its anti-particle
      CT = 7;
      KATT3 = -KATT0;
      KATT2 = KATT3;
      //KATT0=KATT0;
      if (KATT00 > 0) { //qqbar->qqbar
        max_color++;
        color0 = max_color;
        anti_color0 = 0;
        color2 = 0;
        anti_color2 = max_color;
        color3 = 0;
        anti_color3 = backup_color0;
      } else { //qbarq->qbarq
        max_color++;
        color0 = 0;
        anti_color0 = max_color;
        color2 = max_color;
        anti_color2 = 0;
        color3 = backup_anti_color0;
        anti_color3 = 0;
      }
    }
  }
}


void Matter::colljet22(int CT, double temp, double qhat0ud, double v0[4],
                       double p0[4], double p2[4], double p3[4], double p4[4],
                       double &qt) {
  //
  //    p0 initial jet momentum, output to final momentum
  //    p2 final thermal momentum,p3 initial termal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************
  p4[1] = p0[1];
  p4[2] = p0[2];
  p4[3] = p0[3];
  p4[0] = p0[0];
  //************************************************************

  //    transform to local comoving frame of the fluid
  //  cout << endl;
  //  cout << "flow  "<< v0[1] << " " << v0[2] << " " << v0[3] << " "<<" Elab " << p0[0] << endl;

  trans(v0, p0);
  //  cout << p0[0] << " " << sqrt(qhat0ud) << endl;

  //  cout << sqrt(pow(p0[1],2)+pow(p0[2],2)+pow(p0[3],2)) << " " << p0[1] << " " << p0[2] << " " << p0[3] << endl;

  //************************************************************
  trans(v0, p4);
  //************************************************************

  //    sample the medium parton thermal momentum in the comoving frame

  double xw;
  double razim;
  double rcos;
  double rsin;

  double ss;
  double tmin;
  double tmid;
  double tmax;

  double rant;
  double tt;

  double uu;
  double ff = 0.0;
  double rank;

  double mmax;
  double msq = 0.0;

  double f1;
  double f2;

  double p0ex[4] = {0.0};
  double vc[4] = {0.0};

  int ct1_loop, ct2_loop, flag1, flag2;

  flag1 = 0;
  flag2 = 0;

  //  Initial 4-momentum of jet
  //
  //************************************************************
  p4[1] = p0[1];
  p4[2] = p0[2];
  p4[3] = p0[3];
  p4[0] = p0[0];
  //************************************************************

  int ic = 0;

  ct1_loop = 0;
  do {
    ct1_loop++;
    if(flag2 == 1 || ct1_loop > 1e6){
       flag1 = 1;
       break;
    }
    ct2_loop = 0;
    do {
      ct2_loop++;
      if(ct2_loop > 1e6){
         flag2 = 1;
         break;
      }
      xw = 15.0 * ran0(&NUM1);
      razim = 2.0 * pi * ran0(&NUM1);
      rcos = 1.0 - 2.0 * ran0(&NUM1);
      rsin = sqrt(1.0 - rcos * rcos);
      //
      p2[0] = xw * temp;
      p2[3] = p2[0] * rcos;
      p2[1] = p2[0] * rsin * cos(razim);
      p2[2] = p2[0] * rsin * sin(razim);

      //
      //    cms energy
      //
      ss =
          2.0 * (p0[0] * p2[0] - p0[1] * p2[1] - p0[2] * p2[2] - p0[3] * p2[3]);

      //	if(ss.lt.2.d0*qhat0ud) goto 14

      tmin = qhat0ud;
      tmid = ss / 2.0;
      tmax = ss - qhat0ud;

      //    use (s^2+u^2)/(t+qhat0ud)^2 as scattering cross section in the
      //
      rant = ran0(&NUM1);
      tt = rant * ss;

      //		ic+=1;
      //		cout << p0[0] << "  " << p2[0] <<  endl;
      //		cout << tt << "  " << ss <<  "" << qhat0ud <<endl;
      //		cout << ic << endl;

    } while ((tt < qhat0ud) || (tt > (ss - qhat0ud)));

    f1 = pow(xw, 3) / (exp(xw) - 1) / 1.4215;
    f2 = pow(xw, 3) / (exp(xw) + 1) / 1.2845;

    uu = ss - tt;

    if (CT == 1) {
      ff = f1;
      mmax =
          4.0 / pow(ss, 2) *
          (3.0 - tmin * (ss - tmin) / pow(ss, 2) +
           (ss - tmin) * ss / pow(tmin, 2) + tmin * ss / pow((ss - tmin), 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (3.0 - tt * uu / pow(ss, 2) + uu * ss / pow(tt, 2) +
             tt * ss / pow(uu, 2)) /
            mmax;
    }

    if (CT == 2) {
      ff = f1;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / tmin /
                  (ss - tmin) -
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / tt / uu -
             (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) /
            (mmax + 4.0);
    }

    if (CT == 3) {
      ff = f2;
      if (((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
           4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin)) >
          ((pow(ss, 2) + pow((ss - tmax), 2) / pow(tmax, 2) +
            4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss /
                (ss - tmax)))) {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin));
      } else {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmax), 2)) / pow(tmax, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss / (ss - tmax));
      }
      //
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / ss / uu) /
            mmax;
    }

    if (CT == 13) {
      ff = f1;

      if (((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
           4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin)) >
          ((pow(ss, 2) + pow((ss - tmax), 2) / pow(tmax, 2) +
            4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss /
                (ss - tmax)))) {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / ss / (ss - tmin));
      } else {
        mmax =
            4.0 / pow(ss, 2) *
            ((pow(ss, 2) + pow((ss - tmax), 2)) / pow(tmax, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmax), 2)) / ss / (ss - tmax));
      }
      //
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
             4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / ss / uu) /
            mmax;
    }

    if (CT == 4) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(ss, 2) + pow(uu, 2)) / pow(tt, 2)) / mmax;
    }

    if (CT == 5) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
              (pow(ss, 2) + pow(tmin, 2)) / pow((ss - tmin), 2) -
              2.0 / 3.0 * pow(ss, 2) / tmin / (ss - tmin));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 *
             ((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2) +
              (pow(ss, 2) + pow(tt, 2)) / pow(uu, 2) -
              2.0 / 3.0 * pow(ss, 2) / tt / uu)) /
            mmax;
    }

    if (CT == 6) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) / (mmax + 0.5);
    }

    if (CT == 7) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(ss, 2) + pow((ss - tmin), 2)) / pow(tmin, 2) +
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2) +
              2.0 / 3.0 * pow((ss - tmin), 2) / ss / tmin);
      msq = (pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
             (4.0 / 9.0 *
              (((pow(ss, 2) + pow(uu, 2)) / pow(tt, 2)) +
               (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2) +
               2.0 / 3.0 * pow(uu, 2) / ss / tt))) /
            mmax;
    }

    if (CT == 8) {
      ff = f2;
      mmax = 4.0 / pow(ss, 2) *
             (4.0 / 9.0 * (pow(tmin, 2) + pow((ss - tmin), 2)) / tmin /
                  (ss - tmin) -
              (pow(tmin, 2) + pow((ss - tmin), 2)) / pow(ss, 2));
      msq = pow((1.0 / p0[0] / p2[0] / 2.0), 2) *
            (4.0 / 9.0 * (pow(tt, 2) + pow(uu, 2)) / tt / uu -
             (pow(tt, 2) + pow(uu, 2)) / pow(ss, 2)) /
            (mmax + 4.0);
    }

    rank = ran0(&NUM1);
  } while (rank > (msq * ff));

  if(flag1 == 1 || flag2 == 1){ // scatterings cannot be properly sampled
    transback(v0, p0);
    transback(v0, p4);
    qt = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 0;
    p2[3] = 0;
    p3[0] = 0;
    p3[1] = 0;
    p3[2] = 0;
    p3[3] = 0;
    return;
  }

  //
  p3[1] = p2[1];
  p3[2] = p2[2];
  p3[3] = p2[3];
  p3[0] = p2[0];

  //    velocity of the center-of-mass
  //
  vc[1] = (p0[1] + p2[1]) / (p0[0] + p2[0]);
  vc[2] = (p0[2] + p2[2]) / (p0[0] + p2[0]);
  vc[3] = (p0[3] + p2[3]) / (p0[0] + p2[0]);
  //
  //    transform into the cms frame
  //
  trans(vc, p0);
  trans(vc, p2);
  //
  //    cm momentum
  //
  double pcm = p2[0];
  //
  //    sample transverse momentum transfer with respect to jet momentum
  //    in cm frame
  //
  double ranp = 2.0 * pi * ran0(&NUM1);
  //
  //    transverse momentum transfer
  //
  qt = sqrt(pow(pcm, 2) - pow((tt / 2.0 / pcm - pcm), 2));
  double qx = qt * cos(ranp);
  double qy = qt * sin(ranp);

  //
  //    longitudinal momentum transfer
  //
  double qpar = tt / 2.0 / pcm;
  //
  //    qt is perpendicular to pcm, need to rotate back to the cm frame
  //
  double upt = sqrt(p2[1] * p2[1] + p2[2] * p2[2]) / p2[0];
  double upx = p2[1] / p2[0];
  double upy = p2[2] / p2[0];
  double upz = p2[3] / p2[0];
  //
  //    momentum after collision in cm frame
  //
  p2[1] = p2[1] - qpar * upx;
  p2[2] = p2[2] - qpar * upy;
  if (upt != 0.0) {
    p2[1] = p2[1] + (upz * upx * qy + upy * qx) / upt;
    p2[2] = p2[2] + (upz * upy * qy - upx * qx) / upt;
  }
  p2[3] = p2[3] - qpar * upz - upt * qy;

  p0[1] = -p2[1];
  p0[2] = -p2[2];
  p0[3] = -p2[3];
  //
  //    transform from cm back to the comoving frame
  //
  transback(vc, p2);
  transback(vc, p0);

  //************************************************************
  //
  //     calculate qt in the rest frame of medium
  //
  //  if(p0[4]>p2[4])
  //	{
  rotate(p4[1], p4[2], p4[3], p0, 1);
  qt = sqrt(pow(p0[1], 2) + pow(p0[2], 2));
  rotate(p4[1], p4[2], p4[3], p0, -1);
  //	}
  //  else
  //	{
  //	  rotate(p4[1],p4[2],p4[3],p2,1);
  //	  qt=sqrt(pow(p2[1],2)+pow(p2[2],2));
  //	  rotate(p4[1],p4[2],p4[3],p2,-1);
  //	}
  //************************************************************

  //
  //    transform from comoving frame to the lab frame
  //
  transback(v0, p2);
  transback(v0, p0);
  transback(v0, p3);

  //************************************************************
  transback(v0, p4);
  //************************************************************
}

//.........................................................................
void Matter::collHQ22(int CT, double temp, double qhat0ud, double v0[4],
                      double p0[4], double p2[4], double p3[4], double p4[4],
                      double &qt) {
  //
  //    HQ 2->2 scatterings
  //    p0 initial HQ momentum, output to final momentum
  //    p2 final thermal momentum, p3 initial thermal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************

  // transform to local comoving frame of the fluid
  trans(v0, p0);

  //************************************************************

  //    sample the medium parton thermal momentum in the comoving frame

  double xw;
  double razim;
  double rcos;
  double rsin;

  double ss;

  double rant;
  double tt;

  double uu;
  double ff = 0.0;
  double rank;

  double msq = 0.0;

  double e2, theta2, theta4, phi24; // the four independent variables
  double e1, e4, p1, cosTheta24, downFactor,
      sigFactor; // other useful variables
  double HQmass, fBmax, fFmax, fB, fF, maxValue;
  int index_p1, index_T, index_e2;
  int ct1_loop, ct2_loop, flag1, flag2;

  flag1 = 0;
  flag2 = 0;

  // continue this function for HQ scattering

  HQmass = p0[0] * p0[0] - p0[1] * p0[1] - p0[2] * p0[2] - p0[3] * p0[3];
  if (HQmass > 1e-12) {
    HQmass = sqrt(HQmass);
  } else {
    HQmass = 0.0;
  }

  //    Initial 4-momentum of HQ
  //
  //************************************************************
  p4[1] = p0[1];
  p4[2] = p0[2];
  p4[3] = p0[3];
  p4[0] = p0[0];
  //************************************************************

  p1 = sqrt(p0[1] * p0[1] + p0[2] * p0[2] + p0[3] * p0[3]);
  index_p1 = (int)((p1 - min_p1) / bin_p1);
  index_T = (int)((temp - min_T) / bin_T);
  if (index_p1 >= N_p1) {
    index_p1 = N_p1 - 1;
    cout << "warning: p1 is over p_max: " << p1 << endl;
  }
  if (index_T >= N_T) {
    index_T = N_T - 1;
    cout << "warning: T is over T_max: " << temp << endl;
  }
  if (index_T < 0) {
    index_T = 0;
    cout << "warning: T is below T_min: " << temp << endl;
  }

  fBmax = distFncBM[index_T][index_p1];
  fFmax = distFncFM[index_T][index_p1]; // maximum of f(xw) at given p1 and T

  maxValue = 10.0; // need actual value later

  ct1_loop = 0;
  do { // sample p2 (light parton) using distribution integrated over 3 angles
    ct1_loop++;
    if (ct1_loop > 1e6) {
      //            cout << "cannot sample light parton for HQ scattering ..." << endl;
      flag1 = 1;
      break;
    }
    xw = max_e2 * ran0(&NUM1);
    index_e2 = (int)((xw - min_e2) / bin_e2);
    if (index_e2 >= N_e2)
      index_e2 = N_e2 - 1;
    if (CT == 11) { // qc->qc
      ff = distFncF[index_T][index_p1][index_e2] / fFmax;
      maxValue = distMaxF[index_T][index_p1][index_e2];
    } else if (CT == 12) { // gc->gc
      ff = distFncB[index_T][index_p1][index_e2] / fBmax;
      maxValue = distMaxB[index_T][index_p1][index_e2];
    } else {
      cout << "Wrong HQ channel ID" << endl;
      exit(EXIT_FAILURE);
    }
  } while (ran0(&NUM1) > ff);

  e2 = xw * temp;
  e1 = p0[0];

  // now e2 is fixed, need to sample the remaining 3 variables
  ct2_loop = 0;
  do {
    ct2_loop++;
    if (ct2_loop > 1e6) {
      cout << "cannot sample final states for HQ scattering ..." << endl;
      flag2 = 1;
      break;
    }

    theta2 = pi * ran0(&NUM1);
    theta4 = pi * ran0(&NUM1);
    phi24 = 2.0 * pi * ran0(&NUM1);

    cosTheta24 =
        sin(theta2) * sin(theta4) * cos(phi24) + cos(theta2) * cos(theta4);
    downFactor = e1 - p1 * cos(theta4) + e2 - e2 * cosTheta24;
    e4 = (e1 * e2 - p1 * e2 * cos(theta2)) / downFactor;
    sigFactor = sin(theta2) * sin(theta4) * e2 * e4 / downFactor;

    // calculate s,t,u, different definition from light quark -- tt, uu are negative
    ss = 2.0 * e1 * e2 + HQmass * HQmass - 2.0 * p1 * e2 * cos(theta2);
    tt = -2.0 * e2 * e4 * (1.0 - cosTheta24);
    uu = 2.0 * HQmass * HQmass - ss - tt;

    // re-sample if the kinematic cuts are not satisfied
    if (ss <= 2.0 * qhat0ud || tt >= -qhat0ud || uu >= -qhat0ud) {
      rank = ran0(&NUM1);
      sigFactor = 0.0;
      msq = 0.0;
      continue;
    }

    if (CT == 11) { // qc->qc
      ff =
          (1.0 / (exp(e2 / temp) + 1.0)) * (1.0 - 1.0 / (exp(e4 / temp) + 1.0));
      sigFactor = sigFactor * ff;
      msq = Mqc2qc(ss, tt, HQmass) / maxValue;
    }

    if (CT == 12) { // gc->gc
      ff =
          (1.0 / (exp(e2 / temp) - 1.0)) * (1.0 + 1.0 / (exp(e4 / temp) - 1.0));
      sigFactor = sigFactor * ff;
      msq = Mgc2gc(ss, tt, HQmass) / maxValue;
    }

    rank = ran0(&NUM1);

  } while (rank > (msq * sigFactor));

  if (flag1 == 0 && flag2 == 0) {

    // pass p2 value to p3 for initial thermal parton
    p3[1] = e2 * sin(theta2);
    p3[2] = 0.0;
    p3[3] = e2 * cos(theta2);
    p3[0] = e2;

    // calculate momenta of outgoing particles
    // here p2 is for p4 (light parton) in my note

    p2[1] = e4 * sin(theta4) * cos(phi24);
    p2[2] = e4 * sin(theta4) * sin(phi24);
    p2[3] = e4 * cos(theta4);
    p2[0] = e4;

    // rotate randomly in xy plane (jet is in z), because p3 is assigned in xz plane with bias
    double th_rotate = 2.0 * pi * ran0(&NUM1);
    double p3x_rotate = p3[1] * cos(th_rotate) - p3[2] * sin(th_rotate);
    double p3y_rotate = p3[1] * sin(th_rotate) + p3[2] * cos(th_rotate);
    double p2x_rotate = p2[1] * cos(th_rotate) - p2[2] * sin(th_rotate);
    double p2y_rotate = p2[1] * sin(th_rotate) + p2[2] * cos(th_rotate);
    p3[1] = p3x_rotate;
    p3[2] = p3y_rotate;
    p2[1] = p2x_rotate;
    p2[2] = p2y_rotate;

    // Because we treated p0 (p1 in my note for heavy quark) as the z-direction, proper rotations are necessary here
    rotate(p4[1], p4[2], p4[3], p2, -1);
    rotate(p4[1], p4[2], p4[3], p3, -1);

    p0[1] = p4[1] + p3[1] - p2[1];
    p0[2] = p4[2] + p3[2] - p2[2];
    p0[3] = p4[3] + p3[3] - p2[3];
    p0[0] =
        sqrt(p0[1] * p0[1] + p0[2] * p0[2] + p0[3] * p0[3] + HQmass * HQmass);

    // Debug
    if (fabs(p0[0] + p2[0] - p3[0] - p4[0]) > 0.00001) {
      cout << "Violation of energy conservation in HQ 2->2 scattering:  "
           << fabs(p0[0] + p2[0] - p3[0] - p4[0]) << endl;
    }

    // calculate qt in the rest frame of medium
    rotate(p4[1], p4[2], p4[3], p0, 1);
    qt = sqrt(pow(p0[1], 2) + pow(p0[2], 2));
    rotate(p4[1], p4[2], p4[3], p0, -1);

    // transform from comoving frame to the lab frame
    transback(v0, p2);
    transback(v0, p0);
    transback(v0, p3);
    transback(v0, p4);

  } else { // no scattering
    transback(v0, p0);
    transback(v0, p4);
    qt = 0;
    p2[0] = 0;
    p2[1] = 0;
    p2[2] = 0;
    p2[3] = 0;
    p3[0] = 0;
    p3[1] = 0;
    p3[2] = 0;
    p3[3] = 0;
  }
}

double Matter::Mqc2qc(double s, double t, double M) {

  double m2m = M * M;
  double u = 2.0 * m2m - s - t;
  double MM;

  MM = 64.0 / 9.0 * (pow((m2m - u), 2) + pow((s - m2m), 2) + 2.0 * m2m * t) /
       t / t;

  return (MM);
}

double Matter::Mgc2gc(double s, double t, double M) {

  double m2m = M * M;
  double u = 2.0 * m2m - s - t;
  double MM;

  MM = 32.0 * (s - m2m) * (m2m - u) / t / t;
  MM = MM + 64.0 / 9.0 * ((s - m2m) * (m2m - u) + 2.0 * m2m * (s + m2m)) /
                pow((s - m2m), 2);
  MM = MM + 64.0 / 9.0 * ((s - m2m) * (m2m - u) + 2.0 * m2m * (u + m2m)) /
                pow((u - m2m), 2);
  MM = MM + 16.0 / 9.0 * m2m * (4.0 * m2m - t) / ((s - m2m) * (m2m - u));
  MM = MM + 16.0 * ((s - m2m) * (m2m - u) + m2m * (s - u)) / (t * (s - m2m));
  MM = MM + 16.0 * ((s - m2m) * (m2m - u) - m2m * (s - u)) / (t * (u - m2m));

  return (MM);
}

void Matter::trans(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] - ppar);
  p[1] = p[1] + v[1] * gavv;
  p[2] = p[2] + v[2] * gavv;
  p[3] = p[3] + v[3] * gavv;
}

void Matter::transback(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (-ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] + ppar);
  p[1] = p[1] - v[1] * gavv;
  p[2] = p[2] - v[2] * gavv;
  p[3] = p[3] - v[3] * gavv;
}

void Matter::rotate(double px, double py, double pz, double pr[4], int icc) {
  //     input:  (px,py,pz), (wx,wy,wz), argument (i)
  //     output: new (wx,wy,wz)
  //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
  //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

  double wx, wy, wz, E, pt, w, cosa, sina, cosb, sinb;
  double wx1, wy1, wz1;

  wx = pr[1];
  wy = pr[2];
  wz = pr[3];

  E = sqrt(px * px + py * py + pz * pz);
  pt = sqrt(px * px + py * py);

  w = sqrt(wx * wx + wy * wy + wz * wz);

  if (pt == 0) {
    cosa = 1;
    sina = 0;
  } else {
    cosa = px / pt;
    sina = py / pt;
  }

  cosb = pz / E;
  sinb = pt / E;

  if (icc == 1) {
    wx1 = wx * cosb * cosa + wy * cosb * sina - wz * sinb;
    wy1 = -wx * sina + wy * cosa;
    wz1 = wx * sinb * cosa + wy * sinb * sina + wz * cosb;
  }

  else {
    wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
    wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
    wz1 = -wx * sinb + wz * cosb;
  }

  wx = wx1;
  wy = wy1;
  wz = wz1;

  pr[1] = wx;
  pr[2] = wy;
  pr[3] = wz;

  //  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);
}

float Matter::ran0(long *idum)

{
  const int IM1 = 2147483563;
  const int IM2 = 2147483399;
  const double AM = (1.0 / IM1);
  const int IMM1 = (IM1 - 1);
  const int IA1 = 40014;
  const int IA2 = 40692;
  const int IQ1 = 53668;
  const int IQ2 = 52774;
  const int IR1 = 12211;
  const int IR2 = 3791;
  const int NTAB = 32;
  const int NDIV = (1 + IMM1 / NTAB);
  const double EPS = 1.2e-7;
  const double RNMX = (1.0 - EPS);

  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1)
      *idum = 1;
    else
      *idum = -(*idum);
    for (j = NTAB + 7; j >= 0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}

//////////////////////////////////////////////////////////////////////////////////////

double Matter::solve_alphas(double var_qhat, double var_ener, double var_temp) {

  double preFactor = 42.0 * Ca * zeta3 / pi;

  // reference: qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*2.0*pi*tempLoc*tempLoc/4.0/muD2);
  double max_qhat =
      preFactor * pow(0.5, 2) * pow(var_temp, 3) *
      log(5.7 * max(var_ener, 2.0 * pi * var_temp) / 24 / pi / 0.5 / var_temp);

  if (max_qhat < var_qhat) {
    JSINFO << "qhat exceeds HTL calculation, use alpha_s = 0.5";
    return (0.5);
  }

  double solution = sqrt(
      var_qhat / preFactor / pow(var_temp, 3) /
      log(5.7 * max(var_ener, 2.0 * pi * var_temp) / 24 / pi / 0.2 / var_temp));
  double fnc_value, fnc_derivative;
  fnc_value = fnc0_alphas(solution, var_qhat, var_ener, var_temp);
  fnc_derivative =
      fnc0_derivative_alphas(solution, var_qhat, var_ener, var_temp);

  //cout << "initial guess: " << solution << "  " << fnc_value << endl;

  while (fabs(fnc_value / var_qhat) > 0.001) {

    solution = solution - fnc_value / fnc_derivative;
    fnc_value = fnc0_alphas(solution, var_qhat, var_ener, var_temp);
    fnc_derivative =
        fnc0_derivative_alphas(solution, var_qhat, var_ener, var_temp);
  }

  if (solution < 0.0 || solution > 0.5) {
    JSINFO << "unreasonable alpha_s: " << solution << " use alpha_s = 0.5";
    solution = 0.5;
  }

  return (solution);
}

double Matter::fnc0_alphas(double var_alphas, double var_qhat, double var_ener,
                           double var_temp) {

  double preFactor = 42.0 * Ca * zeta3 / pi;
  return (preFactor * var_alphas * var_alphas * pow(var_temp, 3) *
              log(5.7 * max(var_ener, 2.0 * pi * var_temp) / 24 / pi /
                  var_alphas / var_temp) -
          var_qhat);
}

double Matter::fnc0_derivative_alphas(double var_alphas, double var_qhat,
                                      double var_ener, double var_temp) {

  double preFactor = 42.0 * Ca * zeta3 / pi;
  return (preFactor * pow(var_temp, 3) *
          (2.0 * var_alphas *
               log(5.7 * max(var_ener, 2.0 * pi * var_temp) / 24 / pi /
                   var_alphas / var_temp) -
           var_alphas));
}

void Matter::read_tables() { // intialize various tables for LBT

  //...read scattering rate
  int it, ie;
  int n = 450;
  //ifstream f1("LBT-tables/ratedata");
  //if(!f1.is_open())
  //  {
  //    cout<<"Erro openning date file1!\n";
  //  }
  //else
  //  {
  //    for(int i=1;i<=n;i++)
  //      {
  //        f1>>it>>ie;
  //        f1>>qhatG[it][ie]>>Rg[it][ie]>>Rg1[it][ie]>>Rg2[it][ie]>>Rg3[it][ie]>>qhatLQ[it][ie]>>Rq[it][ie]>>Rq3[it][ie]>>Rq4[it][ie]>>Rq5[it][ie]>>Rq6[it][ie]>>Rq7[it][ie]>>Rq8[it][ie];
  //      }
  //  }
  //f1.close();

  // duplicate for heavy quark
  ifstream f11("LBT-tables/ratedata-HQ");
  if (!f11.is_open()) {
    cout << "Erro openning HQ data file!\n";
  } else {
    for (int i = 1; i <= n; i++) {
      f11 >> it >> ie;
      f11 >> RHQ[it][ie] >> RHQ11[it][ie] >> RHQ12[it][ie] >> qhatHQ[it][ie];
    }
  }
  f11.close();

  // preparation for HQ 2->2
  ifstream fileB("LBT-tables/distB.dat");
  if (!fileB.is_open()) {
    cout << "Erro openning data file distB.dat!" << endl;
  } else {
    for (int i = 0; i < N_T; i++) {
      for (int j = 0; j < N_p1; j++) {
        double dummy_T, dummy_p1;
        fileB >> dummy_T >> dummy_p1;
        if (fabs(min_T + (0.5 + i) * bin_T - dummy_T) > 1.0e-5 ||
            fabs(min_p1 + (0.5 + j) * bin_p1 - dummy_p1) > 1.0e-5) {
          cout << "Erro in reading data file distB.dat!" << endl;
          exit(EXIT_FAILURE);
        }
        fileB >> distFncBM[i][j];
        for (int k = 0; k < N_e2; k++)
          fileB >> distFncB[i][j][k];
        for (int k = 0; k < N_e2; k++)
          fileB >> distMaxB[i][j][k];
      }
    }
  }
  fileB.close();

  ifstream fileF("LBT-tables/distF.dat");
  if (!fileF.is_open()) {
    cout << "Erro openning data file distF.dat!" << endl;
  } else {
    for (int i = 0; i < N_T; i++) {
      for (int j = 0; j < N_p1; j++) {
        double dummy_T, dummy_p1;
        fileF >> dummy_T >> dummy_p1;
        if (fabs(min_T + (0.5 + i) * bin_T - dummy_T) > 1.0e-5 ||
            fabs(min_p1 + (0.5 + j) * bin_p1 - dummy_p1) > 1.0e-5) {
          cout << "Erro in reading data file distF.dat!" << endl;
          exit(EXIT_FAILURE);
        }
        fileF >> distFncFM[i][j];
        for (int k = 0; k < N_e2; k++)
          fileF >> distFncF[i][j][k];
        for (int k = 0; k < N_e2; k++)
          fileF >> distMaxF[i][j][k];
      }
    }
  }
  fileF.close();
}





