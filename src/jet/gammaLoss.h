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

  //SC: for interface with hydro
  double fillQhatTab(double y);
  double fncQhat(double zeta);
  double fncAvrQhat(double zeta, double tau);

  bool gammaLoss_on, in_vac, brick_med, recoil_on, broadening_on;

  double initR0, initRx, initRy, initRz, initVx, initVy, initVz, initRdotV,
      initVdotV, initEner;

  double tStart;// = 0.6;
  double hydro_Tc, brick_length;
  int iEvent;
  bool debug_flag = 0;
  long NUM1;


  //qhat stuff
  double absFactor(TLorentzVector pVec, double T);
  double absFactor2(TLorentzVector pVec, double T);
  bool isAbsorbed(TLorentzVector pVec, double T, double delTime);

  // flag to make sure initialize only once
  static bool flag_init;

private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<gammaLoss> reg;
};

#endif // gammaLoss_H
