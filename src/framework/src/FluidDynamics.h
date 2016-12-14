// Framework test (dummy) FluidDynamics class (to be changed with real implemenation)

#ifndef FLUIDDYNAMICS_H
#define FLUIDDYNAMICS_H

#include "JetScapeModuleBase.h"

class FluidDynamics : public JetScapeModuleBase //, std::enable_shared_from_this<FluidDynamics>
{
  
 public:
  
  FluidDynamics();
  FluidDynamics(string m_name) : JetScapeModuleBase (m_name)
    {eta=-99.99; SetId("FluidDynamics");}
  virtual ~FluidDynamics();

  virtual void Init();
  virtual void Exec();
  
  void SetEtat(double m_eta) {eta=m_eta;}
  double GetQhat() {return eta;}

  // slots for "jet" signals
  void UpdateEnergyDeposit(int t, double edop);
  void GetEnergyDensity(int t,double& edensity);
  
 private:

  double eta;
};

#endif
