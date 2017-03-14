//  SignalManager instance class (meant as singelton)

#ifndef JETSCAPESIGNALMANAGER_H
#define JETSCAPESIGNALMANAGER_H

#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "HardProcess.h"
#include "JetScapeWriter.h"

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
  weak_ptr<FluidDynamics> GetHydroPointer() {return hydro;}
  
  void SetJetEnergyLossManagerPointer(shared_ptr<JetEnergyLossManager> m_jloss) {jloss=m_jloss;}
  weak_ptr<JetEnergyLossManager> GetJetEnergyLossManagerPointer() {return jloss;}
  
  void SetHardProcessPointer(shared_ptr<HardProcess> m_hardp) {hardp=m_hardp;}
  weak_ptr<HardProcess> GetHardProcessPointer() {return hardp;}
  
  void SetWriterPointer(shared_ptr<JetScapeWriter> m_writer) {writer=m_writer;}
  weak_ptr<JetScapeWriter> GetWriterPointer() {return writer;}
  
  void ConnectJetSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j);
  void ConnectAddJetSourceSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetTemperatureSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetHydroCellSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetHardPartonListSignal(shared_ptr<JetEnergyLossManager> jm);
    
  void DisconnectSignal() {}; // to be implememted if needed ...

  void CleanUp();
  
  int GetNumberOfJetSignals() {return num_jet_signals;}
  int GetNumberOfEdensitySignals() {return num_edensity_signals;}
  int GetNumberOfAddJetSourceSignals() {return num_AddJetSourceSignals;}
  int GetNumberOfGetTemperatureSignals() {return num_GetTemperatureSignals;}
  int GetNumberOfGetHydroCellSignals() {return num_GetHydroCellSignals;}
  
  void PrintJetSignalMap();
  void PrintEdensitySignalMap();
  void PrintAddJetSourceSignalMap() {};
  void PrintGetTemperatureSignalMap() {};
  void PrintGetHydroCellSignalMap();
  
 private:

  JetScapeSignalManager() {};
  JetScapeSignalManager(JetScapeSignalManager const&) {};
  static JetScapeSignalManager* m_pInstance;

  weak_ptr<FluidDynamics> hydro;
  weak_ptr<JetEnergyLossManager> jloss;
  weak_ptr<HardProcess> hardp;
  weak_ptr<JetScapeWriter> writer;
  
  int num_jet_signals=0;
  int num_edensity_signals=0;
  int num_AddJetSourceSignals=0;
  int num_GetTemperatureSignals=0;
  int num_GetHydroCellSignals=0;

  map<int,weak_ptr<JetEnergyLoss>> jet_signal_map;
  map<int,weak_ptr<JetEnergyLoss>> edensity_signal_map;
  map<int,weak_ptr<JetEnergyLoss>> AddJetSourceSignal_map;
  map<int,weak_ptr<JetEnergyLoss>> GetTemperatureSignal_map;
  map<int,weak_ptr<JetEnergyLoss>> GetHydroCellSignal_map;
  
};

#endif

