// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetEnergyLossManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include <string>
#include "ELossModulesTest.h"
#include <iostream>

using namespace std;

JetEnergyLossManager::JetEnergyLossManager()
{
  SetId("JLossManager");
  GetHardPartonListConnected=false;
  VERBOSE(8);
}

JetEnergyLossManager::~JetEnergyLossManager()
{
  // Check if this is all really needed with shared_ptr ...
  DEBUG;
  Clear();
 
  if (GetNumberOfTasks()>0)
    EraseTaskLast();
}

void JetEnergyLossManager::Clear()
{
  DEBUG<<"Hard Parton List ...";
  
  hp.clear();
  
  int n=GetNumberOfTasks();
  for (int i=1;i<n;i++)
    EraseTaskLast();

  // Clean Up not really working with iterators (see also above!!!) Some logic not clear for me.
  JetScapeSignalManager::Instance()->CleanUp();
  JetScapeTask::ClearTasks();
  
  VERBOSE(8)<<hp.size();
}
  

void JetEnergyLossManager::Init()
{
  INFO<<"Intialize JetEnergyLoss Manager ...";

  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss Manager modules found ...";
      exit(-1);
    }
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  INFO<<"Connect JetEnergyLossManager Signal to Hard Process ...";
  JetScapeSignalManager::Instance()->ConnectGetHardPartonListSignal(shared_from_this());
}


void JetEnergyLossManager::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLossManager::Exec()
{
  INFO<<"Run JetEnergyLoss Manager ...";

  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss Manager modules found ...";
      exit(-1);
    }

  if (GetGetHardPartonListConnected())
    {
      GetHardPartonList(hp);
      
      INFO<<"Number of Hard Partons = "<<hp.size();
      
      for (int i=1;i<hp.size();i++)
	{
	  DEBUG<<"Create the "<<i<<" th copy because number of intital hard partons = "<<hp.size();
	 
	  Add(make_shared<JetEnergyLoss>(*dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0))));
	}
    }
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  DEBUG<<"Check and Create Signal/Slots via JetScapeSignalManaher instance if needed ...";
  
  CreateSignalSlots();

  if (GetGetHardPartonListConnected())
    {
      int n=0;
      for (auto it : GetTaskList())
	{
	  dynamic_pointer_cast<JetEnergyLoss>(it)->AddShowerInitiatingParton(hp[n]);
	  n++;
	}
    }
  
  JetScapeTask::ExecuteTasks();
}

void JetEnergyLossManager::CreateSignalSlots()
{
  for (auto it : GetTaskList())
    for (auto it2 : it->GetTaskList())
      {
	if (!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetJetSignalConnected())
	  JetScapeSignalManager::Instance()->ConnectJetSignal(dynamic_pointer_cast<JetEnergyLoss>(it2));
	
	if (!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetEdensitySignalConnected())	  
	  JetScapeSignalManager::Instance()->ConnectEdensitySignal(dynamic_pointer_cast<JetEnergyLoss>(it2));
	
	if (!dynamic_pointer_cast<JetEnergyLoss>(it2)-> GetGetHydroCellSignalConnected())	  
	  JetScapeSignalManager::Instance()->ConnectGetHydroCellSignal(dynamic_pointer_cast<JetEnergyLoss>(it2));

	// between eloss modules and eloss
	// check the signals itself, probably best via manager in the long run ...
	if(!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetSentInPartonsConnected())
	  JetScapeSignalManager::Instance()->ConnectSentInPartonsSignal(dynamic_pointer_cast<JetEnergyLoss>(it),dynamic_pointer_cast<JetEnergyLoss>(it2));
      }

  JetScapeSignalManager::Instance()->PrintGetHydroCellSignalMap();
  VERBOSE(8);
  JetScapeSignalManager::Instance()->PrintSentInPartonsSignalMap();
}
