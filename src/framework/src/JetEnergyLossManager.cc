// Framework test JetEnergyLossManager class implementation

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
  VERBOSE(8);
}

JetEnergyLossManager::~JetEnergyLossManager()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Destructor called."<<endl;
  VERBOSE(8);
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
}

void JetEnergyLossManager::Exec()
{
  INFO<<"Run JetEnergyLoss Manager ...";

  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss Manager modules found ...";
      exit(-1);
    }

  // --------------------
  // For test only ....
  // mimick at most di-jets here ...
  // add function call to pythia to determine how many partons/jets ...
  if (GetNumberOfTasks()<2)
    {
      INFO<<"Assume Di-Jet : Copy original Eloss Module and add to Manager list ...";
      
      auto jloss=dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0));
      auto jloss2=make_shared<JetEnergyLoss> (*jloss);
      
      Add(jloss2);
    }

  //DEBUG<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  //ClearTaskList();
  //DEBUG<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
   
  //cout<<GetTaskList().size()<<endl;
   // --------------------
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  //cout<<GetNumSignals()<<endl;
  DEBUG<<"Check and Create Signal/Slots via JetScapeSignalManaher instance if needed ...";
  CreateSignalSlots();

  /*
   if (GetNumberOfTasks()>2 && GetNumberOfTasks()<4)
    {
      DEBUG<<"Erase last Task ...";
      EraseTaskLast();
    }

   // Only needed for bookkeeping cleanly and if something gets deleted after connections
   // so if we want to manage only 1 ++ (and delete everyting > 1 then this function is needed!!!)
   JetScapeSignalManager::Instance()->CleanUp();
  */
  
  JetScapeTask::ExecuteTasks();
}

/*
int JetEnergyLossManager::GetNumSignals()
{
  int num_sig=0;
  if (GetNumberOfTasks()>0)
    num_sig=GetNumberOfTasks()*(GetTaskAt(0)->GetNumberOfTasks());
  return num_sig;
}
*/

void JetEnergyLossManager::CreateSignalSlots()
{
  //int i=0;
  for (auto it : GetTaskList())
    {
      for (auto it2 : it->GetTaskList())
	{
	  //cout<<i<<" "<<dynamic_pointer_cast<JetEnergyLoss>(it2)->GetJetSignalConnected()<<endl;
	  //i++;
	  if (!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetJetSignalConnected())
	    {
	      //dynamic_pointer_cast<JetEnergyLoss>(it2)->jetSignal.connect(GetHydroPointer().get(),&FluidDynamics::UpdateEnergyDeposit);
	      //dynamic_pointer_cast<JetEnergyLoss>(it2)->SetJetSignalConnected(true);
	      JetScapeSignalManager::Instance()->ConnectJetSignal(dynamic_pointer_cast<JetEnergyLoss>(it2));
	    }

	  if (!dynamic_pointer_cast<JetEnergyLoss>(it2)->GetEdensitySignalConnected())
	    {
	      JetScapeSignalManager::Instance()->ConnectEdensitySignal(dynamic_pointer_cast<JetEnergyLoss>(it2));
	    }
	}
    }
}
