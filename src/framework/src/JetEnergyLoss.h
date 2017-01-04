// Framework test (dummy) JetEnergyLoss class (to be changed with real implemenation)

#ifndef JETENERGYLOSS_H
#define JETENERGYLOSS_H

#include "JetScapeModuleBase.h"

class JetEnergyLoss : public JetScapeModuleBase //, std::enable_shared_from_this<JetEnergyLoss>
{
  
 public:
  
  JetEnergyLoss();
  JetEnergyLoss(string m_name) : JetScapeModuleBase (m_name)
  {qhat=-99.99;SetId("JetEnergyLoss");}
  JetEnergyLoss(const JetEnergyLoss &j);
  virtual ~JetEnergyLoss();

  virtual void Init();
  virtual void Exec();

  // test only ...
  sigslot::signal2<int, double,multi_threaded_local> jetSignal;
  sigslot::signal2<int, double&,multi_threaded_local> edensitySignal;
  
  void SetQhat(double m_qhat) {qhat=m_qhat;}
  const double GetQhat() const {return qhat;}
  void SetJetSignalConnected(bool m_jetSignalConnected) {jetSignalConnected=m_jetSignalConnected;}
  const bool GetJetSignalConnected() const {return jetSignalConnected;}
  void SetEdensitySignalConnected(bool m_edensitySignalConnected) {edensitySignalConnected=m_edensitySignalConnected;}
  const bool GetEdensitySignalConnected() const {return edensitySignalConnected;}
  
 private:

  double qhat;
  bool jetSignalConnected;
  bool edensitySignalConnected;

};

#endif
