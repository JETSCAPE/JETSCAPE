// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// JetScapeTask class implementation
// Implement truly recursive (not yet done)
// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)

#include "JetScapeTask.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"

#include "JetEnergyLoss.h"

#include <iostream>

using namespace std;

namespace Jetscape {

JetScapeTask::JetScapeTask()
{
  active_exec=true;
  id="";
  my_task_number_ = JetScapeTaskSupport::Instance()->RegisterTask();
  VERBOSE(9);
}

JetScapeTask::~JetScapeTask()
{
  VERBOSE(9);
  DEBUG << "Deleting task with id=" << GetId() << " and TaskNumber= " << get_my_task_number();
}

void JetScapeTask::Init()
{
  DEBUG;
}

void JetScapeTask::InitTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
 
  //In short ...
  for (auto it : tasks)
    it->Init();
}

void JetScapeTask::Exec()
{
  VERBOSE(7);
}

void JetScapeTask::ExecuteTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  if (active_exec)
    for (auto it : tasks)
      it->Exec();
}

void JetScapeTask::ClearTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  for (auto it : tasks)
    it->Clear();
}

void JetScapeTask::WriteTasks(weak_ptr<JetScapeWriter> w)
{
  //VERBOSE(10);
  INFO<<" writer active? " << (w.lock()==NULL ? 0 : 1);
  if (active_exec)
    {
      for (auto it : tasks)
	it->WriteTask(w);
    }
}

void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks)
{
  tasks.push_back(m_tasks);
}

} // end namespace Jetscape
