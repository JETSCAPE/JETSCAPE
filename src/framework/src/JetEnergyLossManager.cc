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
#include <vector>

#include <thread>        
//#include <future>

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
  DEBUG<<"Task Id = "<<this_thread::get_id();
  
  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss Manager modules found ...";
      exit(-1);
    }

  // ----------------------------------
  // Create needed copies and connect signal/slots accordingly ...
  
  if (GetGetHardPartonListConnected())
    {
      GetHardPartonList(hp);
      
      INFO<<" Number of Hard Partons = "<<hp.size();
      
      for (int i=1;i<hp.size();i++)
	{
	  DEBUG<<"Create the "<<i<<" th copy because number of intital hard partons = "<<hp.size();
	 
	  Add(make_shared<JetEnergyLoss>(*dynamic_pointer_cast<JetEnergyLoss>(GetTaskAt(0))));
	}
    }
  
  INFO<<" Found "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules Execute them ... ";
  DEBUG<<"Check and Create Signal/Slots via JetScapeSignalManaher instance if needed ...";
  
  CreateSignalSlots();

  // ----------------------------------
  // Copy Shower initiating partons to JetEnergyLoss tasks ...
  
  if (GetGetHardPartonListConnected())
    {
      int n=0;
      for (auto it : GetTaskList())
	{
	  dynamic_pointer_cast<JetEnergyLoss>(it)->AddShowerInitiatingParton(hp[n]);
	  n++;
	}
    }

  // ----------------------------------
  // quick and dirty here, only include after further testing (flag in init xml files ...)
  bool multiTask=false;

  // ----------------------------------
  //Excute JetEnergyLoss tasks and their subtasks (done via signal/slot) by hand ...
  //needed if only JetEnergyloss tasks in parallel (otherwise could be done via JetScapeTask, then every task a new thread for example)
  // Of course now the number of threads if plainly given by the number of initial hard partons and can exceed the allowed # of threads >1000
  // so to really do this properly one has to think about and limit to number of CPU's * N threads or something in that directions. See a quick attempt below.

  //DEBUG:
  //unsigned num_cpus = thread::hardware_concurrency();
  //cout << "Num of CPU's = " << num_cpus << " threads\n";
  
  if (multiTask)
    {
      int nTasks=GetNumberOfTasks();
      int nCPUs=thread::hardware_concurrency();
      
      vector<thread> threads;
      
      int nMaxThreads=nCPUs*2;
      int n=0;

      INFO<<" Use multi-threading: (max) # of threads = # of CPU's "<<nCPUs<<" (found) * 2";
      
      for (auto it :  GetTaskList())
	{
	  if (it->GetActive())
	    {
	      threads.push_back(thread(&JetEnergyLoss::Exec,dynamic_pointer_cast<JetEnergyLoss>(it)));
	      n++;
	    }
	  if (n==nMaxThreads)
	    {
	      //DEBUGTHREAD<<n;	      
	      for (auto& th : threads) th.join();
	      n=0; threads.clear();

	      //DEBUG:
	      //std::this_thread::sleep_for(std::chrono::milliseconds(5000));
	    }
	}

      if (nTasks<nMaxThreads)
	{      
	  for (auto& th : threads) th.join();      
	  threads.clear();
	}
    }
  // ----------------------------------
  else
    // Standard "serial" execution for the JetEnerguLoss (+submodules) task ...
    JetScapeTask::ExecuteTasks();
  
  //Add acheck if the parton shower was actually created for the Modules ....
  INFO<<" "<<GetNumberOfTasks()<<" Eloss Manager Tasks/Modules finished.";
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
