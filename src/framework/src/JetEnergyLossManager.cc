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
  GetHardPartonListConnected=false;
  VERBOSE(8);
}

JetEnergyLossManager::~JetEnergyLossManager()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Destructor called."<<endl;

  // check if this is all? Last signal not removed in Clear ...
  DEBUG;
  Clear();
  // Check if this is all really needed with shared_ptr ...
  
  if (GetNumberOfTasks()>0)
    EraseTaskLast();
}

void JetEnergyLossManager::Clear()
{
  // copy of vector made due to assignment in HardProcess signal ... (change later with pointers ...)
  // or leave copy in order to finish HardProcess and clear earlier ... !???
  // to be discussed.
  //cout<<hp.size()<<endl;

  DEBUG<<"Hard Parton List ...";
  
  //cout<<hp.size()<<endl;
  hp.clear();
  // erease copied tasks ....
  //cout<<GetNumberOfTasks()<<endl;
  int n=GetNumberOfTasks();
  
  for (int i=1;i<n;i++)
    {
      //cout<<"Loop "<<i<<endl;
      EraseTaskLast();
      //cout<<"After first ..."<<endl;
      //cout<<i<<endl;
    }

  //cout<<GetNumberOfTasks()<<endl;
  //cout<<"out of loop"<<endl;

  // Clean Up not really working with iterators (see also above!!!) Some logic not clear for me.
  JetScapeSignalManager::Instance()->CleanUp();
  
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
      
      // Add Check here ....
      /*
	if (hp.size()<1)
	{
	WARN << " : No valid Energy Loss Manager modules found ...";
	exit(-1);
	}
      */
      
      for (int i=1;i<hp.size();i++)
	{
	  DEBUG<<"Create the "<<i<<" th copy because number of intital hard partons = "<<hp.size();
	  //auto jloss=dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0));
	  //auto jloss2=make_shared<JetEnergyLoss> (*jloss);
	  //Add(jloss2);
	  //in short ..
	  Add(make_shared<JetEnergyLoss>(*dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0))));
	}
    }
  
  // --------------------
  // For test only ....
  // mimick at most di-jets here ...
  // add function call to pythia to determine how many partons/jets ...
  /*
    if (GetNumberOfTasks()<4)
    {
    INFO<<"Assume Di-Jet : Copy original Eloss Module and add to Manager list ...";
    
    auto jloss=dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0));
    auto jloss2=make_shared<JetEnergyLoss> (*jloss);
    Add(jloss2);
    }
  */
  
  //DEBUG<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  //ClearTaskList();
  //DEBUG<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  
  //cout<<GetTaskList().size()<<endl;
  // --------------------
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  //cout<<GetNumSignals()<<endl;
  DEBUG<<"Check and Create Signal/Slots via JetScapeSignalManaher instance if needed ...";
  CreateSignalSlots();

  if (GetGetHardPartonListConnected())
    {
      int n=0;
      for (auto it : GetTaskList())
	{
	  //cout<<hp[n]->pid()<<endl;
	  //cout<<hp[n]<<endl;
	  
	  dynamic_pointer_cast<JetEnergyLoss>(it)->AddShowerInitiatingParton(hp[n]);
	  n++;
	}
    }
  
  /*
   if (GetNumberOfTasks()>2 && GetNumberOfTasks()<4)
    {
      DEBUG<<"Erase last Task ...";
      EraseTaskLast();
      // add in task or overload for modules with signals ...!???
      // proably better to put in task ...
      JetScapeSignalManager::Instance()->CleanUp();
      JetScapeSignalManager::Instance()->PrintGetHydroCellSignalMap();
    }
  */
  
    /*
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

	  if (!dynamic_pointer_cast<JetEnergyLoss>(it2)-> GetGetHydroCellSignalConnected())
	    {
	      JetScapeSignalManager::Instance()->ConnectGetHydroCellSignal(dynamic_pointer_cast<JetEnergyLoss>(it2));
	    }
	}
    }

  JetScapeSignalManager::Instance()->PrintGetHydroCellSignalMap();
}
