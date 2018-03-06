// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETENERGYLOSSMANAGER_H
#define JETENERGYLOSSMANAGER_H

#include "JetScapeTask.h"
#include "JetClass.hpp"
#include "sigslot.h"

#include <vector>

namespace Jetscape {
  /** @class Jet energy loss manager.
   */
class JetEnergyLossManager : public JetScapeTask, public std::enable_shared_from_this<JetEnergyLossManager>
{
  
 public:
  /** Default constructor to create a jet energy loss manager. Sets task ID as "JLossManager". Flag GetHardPartonListConnected is set to false.
   */  
  JetEnergyLossManager();

  /** Destructor for the jet energy loss manager.
   */
  virtual ~JetEnergyLossManager();

  /** It initializes the tasks attached to the jet energy loss manager. It also sends a signal to connect the JetEnergyLoss object to the GetHardPartonList() function of the HardProcess class. It can be overridden by other tasks.
      @sa JetScapeSignalManager to understand the implementation of signal slots philosophy. 
   */  
  virtual void Init();

  /** It reads the Hard Patrons list and calls CreateSignalSlots() function. Then, it executes the energy loss tasks attached with the jet energy loss manager. This function also includes the parallel computing feature. It can be overridden by other tasks.
  */
  virtual void Exec();

  /** It erases the tasks attached with the energy loss manager. It can be overridden by other tasks.
   */
  virtual void Clear();

  /** It writes the output information relevant to the jet energy loss tasks/subtasks into a file. It can be overridden by other tasks.
      @param w A pointer of type JetScapeWriter class.
      @sa JetScapeWriter class for further information. 
  */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
  int GetNumSignals();
  
  /** Uses philosophy of signal slots. Checks whether the attached task is connected via signal slots to the functions UpdateEnergyDeposit(), GetEnergyDensity(), GetHydroCell() (defined in FluidDynamics class), and DoEnergyLoss() (defined in JetEnergyLoss class). If not, then, it sends a signal to these functions.
      @sa JetScapeSignalManager to understand the implementation of signal slots philosophy.
   */
  void CreateSignalSlots();

  /** A signal to connect the JetEnergyLossManager to the function GetHardPartonList() of the class HardProcess. 
   */
  sigslot::signal1<vector<shared_ptr<Parton>>& > GetHardPartonList;

  /** Use the flag m_GetHardPartonListConnected as true, if JetEnergyLossManager had sent a signal to function GetHardPartonList() of the class HardProcess.
      @param m_GetHardPartonListConnected A boolean flag.
   */
  void SetGetHardPartonListConnected(bool m_GetHardPartonListConnected) {GetHardPartonListConnected=m_GetHardPartonListConnected;}

  /** @return GetHardPartonListConnected A boolean flag. Its status indicates whether JetEnergyLossManager had sent a signal to the function GetHardPartonList() of the class HardProcess. 
   */
  const bool GetGetHardPartonListConnected() {return GetHardPartonListConnected;}
  
 private:

  bool GetHardPartonListConnected;
  vector<shared_ptr<Parton>> hp;
  
};

} // end namespace Jetscape

#endif
