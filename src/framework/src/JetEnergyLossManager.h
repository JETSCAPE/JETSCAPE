// Framework test JetEnergyLossManager class 

#ifndef JETENERGYLOSSMANAGER_H
#define JETENERGYLOSSMANAGER_H

#include "JetScapeTask.h"
#include "FluidDynamics.h"
#include "JetClass.hpp"
#include <vector>

class JetEnergyLossManager : public JetScapeTask, public std::enable_shared_from_this<JetEnergyLossManager>
{
  
 public:
  
  JetEnergyLossManager();
  virtual ~JetEnergyLossManager();
  
  virtual void Init();
  virtual void Exec();
  virtual void Clear();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
  int GetNumSignals();
  
  void CreateSignalSlots();
  //void SetHydroPointer(shared_ptr<FluidDynamics> m_hydro) {hydro=m_hydro;}
  //shared_ptr<FluidDynamics> GetHydroPointer() {return hydro;}
  sigslot::signal1<vector<shared_ptr<Parton>>& > GetHardPartonList;

  void SetGetHardPartonListConnected(bool m_GetHardPartonListConnected) {GetHardPartonListConnected=m_GetHardPartonListConnected;}
  const bool GetGetHardPartonListConnected() {return GetHardPartonListConnected;}
  
 private:

  //shared_ptr<FluidDynamics> hydro;
  bool GetHardPartonListConnected;
  vector<shared_ptr<Parton>> hp; // think about memeory here better just a pointer? Or create own copy to be able to delete HardProcess task !????
  
};

#endif
