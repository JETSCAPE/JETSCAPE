// JetScape XML init reader class (meant as singelton)
//  SignalManager instance class

#ifndef JETSCAPESIGNALMANAGER_H
#define JETSCAPESIGNALMANAGER_H

#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"

#include "FluidDynamics.h"

#include<iostream>
#include<string>
#include<map>
#include "sigslot.h"

using namespace sigslot;
using namespace std;

class JetScapeSignalManager //: public sigslot::has_slots<sigslot::multi_threaded_local>
{
  
 public:

  static JetScapeSignalManager* Instance();

  void SetHydroPointer(shared_ptr<FluidDynamics> m_hydro) {hydro=m_hydro;}
  //shared_ptr<FluidDynamics> GetHydroPointer() {return hydro;}
  weak_ptr<FluidDynamics> GetHydroPointer() {return hydro;}
  void SetJetEnergyLossManagerPointer(shared_ptr<JetEnergyLossManager> m_jloss) {jloss=m_jloss;}
  //shared_ptr<JetEnergyLossManager> GetJetEnergyLossManagerPointer() {return jloss;}
  weak_ptr<JetEnergyLossManager> GetJetEnergyLossManagerPointer() {return jloss;}
  
  //void ConnectJetSignal(shared_ptr<JetScapeModuleBase> j);
  void ConnectJetSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j);
  void DisconnectSignal() {};
  //void Clear();

  void CleanUp();
  
  int GetNumberOfJetSignals() {return num_jet_signals;}
  int GetNumberOfEdensitySignals() {return num_edensity_signals;}
  
  void PrintJetSignalMap();
  void PrintEdensitySignalMap();
  
 private:

  JetScapeSignalManager() {};
  JetScapeSignalManager(JetScapeSignalManager const&) {};
  static JetScapeSignalManager* m_pInstance;

  //shared_ptr<FluidDynamics> hydro;
  //shared_ptr<JetEnergyLossManager> jloss;

  weak_ptr<FluidDynamics> hydro;
  weak_ptr<JetEnergyLossManager> jloss;

  int num_jet_signals=0;
  int num_edensity_signals=0;

  map<int,weak_ptr<JetEnergyLoss>> jet_signal_map;
  map<int,weak_ptr<JetEnergyLoss>> edensity_signal_map;
  
};

#endif

