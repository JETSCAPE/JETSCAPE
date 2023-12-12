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

#ifndef gammaLoss_H
#define gammaLoss_H

#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"
#include "TLorentzVector.h"
#include "random"
#include "RtypesCore.h"
#include "TF1.h"
#include "Math/DistSampler.h"

using namespace Jetscape;

class gammaLoss : public JetEnergyLossModule<
                   gammaLoss> //, public std::enable_shared_from_this<gammaLoss>
{
public:
  gammaLoss();
  virtual ~gammaLoss();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn,
                    vector<Parton> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);
  void Dump_pIn_info(int i, vector<Parton> &pIn);
  void doEmission(vector<Parton> &pIn, vector<Parton> &pOut, double deltaT, double time);

  bool gammaLoss_on, in_vac, brick_med, recoil_on, broadening_on, reabsorption;
  int ratesource;
  int emissionOn;
  double integral1, integral2, x0, x1;
  TF1* thermalpdf;
  ROOT::Math::DistSampler* sampler;

  double initR0, initRx, initRy, initRz, initVx, initVy, initVz, initRdotV,
      initVdotV, initEner;

  double tStart;// = 0.6;
  double hydro_Tc, brick_length;
  int iEvent;
  bool debug_flag = 0;
  long NUM1;
  std::mt19937_64 rng_engine;

  //running couplings
  static double alphaS(double temp);
  static double gS(double temp);

  //absorption stuff
  double absFactor1(TLorentzVector pVec, double T);
  double absFactor2(TLorentzVector pVec, double T);
  bool isAbsorbed(TLorentzVector pVec, double T, double delTime);

  //emission
  int photonsProduced(double cell, double temp);
  Parton makeThermalPhoton(double temp, TVector3 vMed, double position[]);
  static double B(double temp);
  static Double_t dRdx(Double_t *x, Double_t *par);
  static Double_t f(Double_t *x, Double_t *par);
  static Double_t g(Double_t *x, Double_t *par);

  std::mt19937_64& getRandomGenerator() {return rng_engine;}
  double genPhi() {return ((float)rand()/RAND_MAX)*2*pi;}
  double genTheta() {return acos((((float)rand()/RAND_MAX)*2)-1);}

  // flag to make sure initialize only once
  static bool flag_init;

private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<gammaLoss> reg;
};

#endif // gammaLoss_H
