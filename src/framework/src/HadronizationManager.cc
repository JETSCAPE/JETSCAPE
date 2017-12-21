
#include "HadronizationManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include <string>
#include "Hadronization.h"

#include <iostream>
#include <vector>

#include <thread>        
#include <future>

using namespace std;

namespace Jetscape {

HadronizationManager::HadronizationManager()
{
  SetId("HadronizationManager");
  GetFinalPartonListConnected = false;
  VERBOSE(8);
}

HadronizationManager::~HadronizationManager()
{   
  Clear();      
  
  // The tasks are hadronization modules
  if (GetNumberOfTasks()>0)
    EraseTaskLast();
}

void HadronizationManager::Clear()
{
  DEBUG<<"Hadronization task List ...";
  
  hd.clear();

  int n=GetNumberOfTasks();
  for (int i=1;i<n;i++)
    EraseTaskLast();
  
  JetScapeSignalManager::Instance()->CleanUp();
  JetScapeTask::ClearTasks();
  VERBOSE(8)<<hd.size();
}

void HadronizationManager::Init()
{
  INFO<<"Intialize Hadronization Manager ...";

  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Hadronization Manager modules found ...";
      exit(-1);
    }
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Hadronization Manager Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  INFO<<"Connect HadronizationManager Signal to Energy Loss ...";
  JetScapeSignalManager::Instance()->ConnectGetFinalPartonListSignal(shared_from_this());

}

void HadronizationManager::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void HadronizationManager::Exec()
{

  VERBOSE(2)<<"Run Hadronization Manager ...";
  DEBUG<<"Task Id = "<<this_thread::get_id();
  
  if (GetNumberOfTasks()<1)
  {
    WARN << " : No valid Hadronization modules found ...";
    exit(-1);
  }

  CreateSignalSlots();

  if (GetGetFinalPartonListConnected())
  {
    GetFinalPartonList(hd);
    VERBOSE(2)<<" There are "<<hd.size() << " partons ready for recombination";
    for (auto it : GetTaskList())
    {
      dynamic_pointer_cast<Hadronization>(it)->AddInPartons(hd);
    }
    JetScapeTask::ExecuteTasks();
  }
  else
  {
    VERBOSE(2)<<" There are no partons ready for recombination" ;
  }

}

void HadronizationManager::CreateSignalSlots()
{
  for (auto it : GetTaskList())
    for (auto it2 : it->GetTaskList())
      if(!dynamic_pointer_cast<Hadronization>(it2)->GetTransformPartonsConnected())
        JetScapeSignalManager::Instance()->ConnectTransformPartonsSignal(dynamic_pointer_cast<Hadronization>(it),dynamic_pointer_cast<Hadronization>(it2));
  
  JetScapeSignalManager::Instance()->PrintTransformPartonsSignalMap();
}











}













