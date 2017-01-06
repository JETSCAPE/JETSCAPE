// Framework test JetEnergyLossManager class 

#ifndef JETENERGYLOSSMANAGER_H
#define JETENERGYLOSSMANAGER_H

#include "JetScapeTask.h"
#include "FluidDynamics.h"

class JetEnergyLossManager : public JetScapeTask
{
  
 public:
  
  JetEnergyLossManager();
  virtual ~JetEnergyLossManager();
  
  virtual void Init();
  virtual void Exec();

  int GetNumSignals();
  
  void CreateSignalSlots();
  //void SetHydroPointer(shared_ptr<FluidDynamics> m_hydro) {hydro=m_hydro;}
  //shared_ptr<FluidDynamics> GetHydroPointer() {return hydro;}
  
 private:

  //shared_ptr<FluidDynamics> hydro;

};

#endif
