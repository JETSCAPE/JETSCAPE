// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETENERGYLOSS_H
#define JETENERGYLOSS_H

#include "JetScapeModuleBase.h"
#include "fluid_dynamics.h"
#include "JetClass.hpp"
#include "JetScapeWriter.h"
#include "PartonShower.h"
#include "helper.h"
#include <vector>
#include <random>

namespace Jetscape {


class JetEnergyLoss : public JetScapeModuleBase, public std::enable_shared_from_this<JetEnergyLoss> 
{
  
 public:
  
  JetEnergyLoss();
  JetEnergyLoss(string m_name) : JetScapeModuleBase (m_name)
  {qhat=-99.99;SetId("JetEnergyLoss");jetSignalConnected=false;
    edensitySignalConnected=false; AddJetSourceSignalConnected=false;
    GetTemperatureSignalConnected=false; GetHydroCellSignalConnected=false;}
  JetEnergyLoss(const JetEnergyLoss &j);
  virtual ~JetEnergyLoss();

  virtual shared_ptr<JetEnergyLoss> Clone() const {return nullptr;}  // const = 0;
  
  virtual void Init();
  virtual void Exec() final; // prevents eloss modules from overwrting and missusing
  virtual void WriteTask(weak_ptr<JetScapeWriter> w); 
  virtual void Clear();
  //virtual void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut) {};
  virtual void DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut) {};
  
  // test only ...
  sigslot::signal2<int, double,multi_threaded_local> jetSignal;
  sigslot::signal2<int, double&,multi_threaded_local> edensitySignal;

  // for brick/gubser test ...
  sigslot::signal5<double, double, double, double, JetSource,multi_threaded_local> AddJetSourceSignal;
  sigslot::signal5<double, double, double, double, double&,multi_threaded_local> GetTemperatureSignal;
  sigslot::signal5<double, double, double, double, FluidCellInfo*,multi_threaded_local> GetHydroCellSignal;

  // signal to all energy loss modules ... get intial list and delta T ... (think more !???)
  // test first ...
  // deltaT , criteria , list
  //sigslot::signal4<double, double, const vector<Parton>&, vector<Parton>&, multi_threaded_local> SentInPartons;
  sigslot::signal4<double, double, vector<Parton>&, vector<Parton>&, multi_threaded_local> SentInPartons;
  sigslot::signal1<vector<Parton>&, multi_threaded_local> GetOutPartons; // probably not needed ... do in SentInPartons with return ...
  
  void SetQhat(double m_qhat) {qhat=m_qhat;}
  const double GetQhat() const {return qhat;}
  
  // old test signals ...
  void SetJetSignalConnected(bool m_jetSignalConnected) {jetSignalConnected=m_jetSignalConnected;}
  const bool GetJetSignalConnected() const {return jetSignalConnected;}  
  void SetEdensitySignalConnected(bool m_edensitySignalConnected) {edensitySignalConnected=m_edensitySignalConnected;}
  const bool GetEdensitySignalConnected() const {return edensitySignalConnected;}

  // signals for jetscape (brick test so far, might be extended)
  void SetAddJetSourceSignalConnected(bool m_AddJetSourceSignalConnected) {AddJetSourceSignalConnected=m_AddJetSourceSignalConnected;}
  const bool GetAddJetSourceSignalConnected() {return AddJetSourceSignalConnected;}
  
  void SetGetTemperatureSignalConnected(bool m_GetTemperatureSignalConnected) {GetTemperatureSignalConnected=m_GetTemperatureSignalConnected;}
  const bool GetGetTemperatureSignalConnected() {return GetTemperatureSignalConnected;}
  
  void SetGetHydroCellSignalConnected(bool m_GetHydroCellSignalConnected) {GetHydroCellSignalConnected=m_GetHydroCellSignalConnected;}
  const bool GetGetHydroCellSignalConnected() {return GetHydroCellSignalConnected;}

  void SetSentInPartonsConnected(bool m_SentInPartonsConnected) {SentInPartonsConnected=m_SentInPartonsConnected;}
  const bool GetSentInPartonsConnected() {return SentInPartonsConnected;}

  void SetGetOutPartonsConnected(bool m_GetOutPartonsConnected) {GetOutPartonsConnected=m_GetOutPartonsConnected;}
  const bool GetGetOutPartonsConnected() {return GetOutPartonsConnected;}
  
  void AddShowerInitiatingParton(shared_ptr<Parton> p) {inP=p;}
  shared_ptr<Parton> GetShowerInitiatingParton() {return inP;}  
  
  void PrintShowerInitiatingParton();

  double GetDeltaT() {return deltaT;}
  double GetMaxT() {return maxT;}
  shared_ptr<PartonShower> GetShower() {return pShower;}
  
 private:

  double deltaT;
  double maxT; // quick fix here ...
  
  double qhat;
  //old test signals 
  bool jetSignalConnected;
  bool edensitySignalConnected;
  
  bool AddJetSourceSignalConnected;
  bool GetTemperatureSignalConnected; //probably not needed if everything via HydroCell
  bool GetHydroCellSignalConnected;
  bool SentInPartonsConnected;
  bool GetOutPartonsConnected;
  
  shared_ptr<Parton> inP;
  //unique_ptr<PartonShower> pShower;
  shared_ptr<PartonShower> pShower;
  
  node vStart;
  node vEnd;

  void DoShower();
  
};

} // end namespace Jetscape

#endif
