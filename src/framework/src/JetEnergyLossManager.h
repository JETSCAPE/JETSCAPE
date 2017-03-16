// Framework test JetEnergyLossManager class 

#ifndef JETENERGYLOSSMANAGER_H
#define JETENERGYLOSSMANAGER_H

#include "JetScapeTask.h"
#include "JetClass.hpp"
#include "sigslot.h"

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

  sigslot::signal1<vector<shared_ptr<Parton>>& > GetHardPartonList;

  void SetGetHardPartonListConnected(bool m_GetHardPartonListConnected) {GetHardPartonListConnected=m_GetHardPartonListConnected;}
  const bool GetGetHardPartonListConnected() {return GetHardPartonListConnected;}
  
 private:

  bool GetHardPartonListConnected;
  vector<shared_ptr<Parton>> hp;
  
};

#endif
