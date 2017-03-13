// JetScapeTask class implementation
// Implement truly recursive (not yet done)
// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)

#include "JetScapeTask.h"
#include "JetScapeLogger.h"
//#include "JetEnergyLoss.h"
//#include "JetEnergyLossManager.h"

#include <iostream>

using namespace std;

JetScapeTask::JetScapeTask()
{
  //Simple Debug replace --> logger
  //cout<<"JetScapeTask : Default Constructor called."<<endl;
  active_exec=true;
  id="";
  VERBOSE(9);
}

JetScapeTask::~JetScapeTask()
{
  //Simple Debug replace --> logger
  //cout<<"JetScapeTask : Default Destructor called."<<endl;
  VERBOSE(9);
}

void JetScapeTask::Init()
{
  DEBUG;
}

void JetScapeTask::InitTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  /*
  for (list<shared_ptr<JetScapeTask>>::iterator it=tasks.begin(); it != tasks.end(); ++it)
    {
      (*it)->Init();
    }
  //Shorter:
  for (auto it=tasks.begin(); it != tasks.end(); ++it)
    {
      (*it)->Init();
    }
  */
  //Even shorter:
  for (auto it : tasks)
    it->Init();
}

void JetScapeTask::Exec()
{
  VERBOSE(7);
  //for (auto it : tasks)
  //it->Exec();
}

void JetScapeTask::ExecuteTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  if (active_exec)
    {
      for (auto it : tasks)
	it->Exec();
    }
}

void JetScapeTask::ClearTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  for (auto it : tasks)
    it->Clear();
}

void JetScapeTask::WriteTasks(weak_ptr<JetScapeWriter> w)
{
  //VERBOSE(7);
  if (active_exec)
    {
      for (auto it : tasks)
	it->WriteTask(w);
    }
}

void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks)
{
  // not sufficent since Matter inherits from JetEnergyLoss
  // so maybe that inheritence is anyhow not necessary !?
  // Just from task/module not JetEnergyloss !??? (Think about)
  //if (dynamic_pointer_cast<JetEnergyLoss>(m_tasks))
  //{
  //  cout<<"Found"<<" "<<m_tasks->GetNumberOfTasks()<<endl;
  //}

  /*
  string m_id=m_tasks->GetId();
    
  if ((int) m_id.find("JetEnergy")>=0)
    {
      VERBOSE(9)<<"Found"<<" "<<m_id;
      DEBUG<<"Create JetScapeEnergyLossManager ...";
      auto jlossmanager = make_shared<JetEnergyLossManager> ();
      cout<<jlossmanager->GetNumberOfTasks()<<endl;
      //cout<<m_tasks->GetNumberOfTasks()<<endl;
	
      //(vector<shared_ptr<JetScapeTask>> (jlossmanager->GetTaskList())).push_back(m_tasks);   
      //jlossmanager->GetTaskList().push_back(dynamic_pointer_cast<JetEnergyLoss>(m_tasks) m_tasks);
      jlossmanager->GetTaskList().push_back(m_tasks);
      jlossmanager->GetTaskList().push_back(m_tasks);
      //cout<<jlossmanager->GetTaskList().size()<<endl;
      cout<<jlossmanager->GetNumberOfTasks()<<endl;
      tasks.push_back(jlossmanager);     
    }
  else
  */
  
  tasks.push_back(m_tasks);
}
