// Framework test (dummy) JetEnergyLoss class (to be changed with real implemenation)

#ifndef JETENERGYLOSS_H
#define JETENERGYLOSS_H

#include "JetScapeModuleBase.h"
#include "fluid_dynamics.h"
#include "JetClass.hpp"
#include "JetScapeWriter.h" // why include needed here !???

//class JetScapeWriter;

class JetEnergyLoss : public JetScapeModuleBase, public std::enable_shared_from_this<JetEnergyLoss> //check memory !?
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
  //virtual JetEnergyLoss *Clone() const = 0;
  
  virtual void Init();
  virtual void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
  //virtual void Clear();
  
  // test only ...
  sigslot::signal2<int, double,multi_threaded_local> jetSignal;
  sigslot::signal2<int, double&,multi_threaded_local> edensitySignal;

  // for brick/gubser test ...
  sigslot::signal5<double, double, double, double, JetSource,multi_threaded_local> AddJetSourceSignal;
  sigslot::signal5<double, double, double, double, double&,multi_threaded_local> GetTemperatureSignal;
  sigslot::signal5<double, double, double, double, FluidCellInfo*,multi_threaded_local> GetHydroCellSignal;
  
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

  void AddShowerInitiatingParton(shared_ptr<Parton> p) {inP=p;}
  //unique_ptr<Parton> GetShowerInitiatingParton() {return inP;}
  shared_ptr<Parton> GetShowerInitiatingParton() {return inP;}
  void PrintShowerInitiatingParton();
  
 private:

  double qhat;
  //old test signals 
  bool jetSignalConnected;
  bool edensitySignalConnected;
  
  bool AddJetSourceSignalConnected;
  bool GetTemperatureSignalConnected; //probably not needed if everything via HydroCell
  bool GetHydroCellSignalConnected;

  //unique_ptr<Parton> inP; //shower/eloss intiated Parton ...
  shared_ptr<Parton> inP;
  
};

#endif
