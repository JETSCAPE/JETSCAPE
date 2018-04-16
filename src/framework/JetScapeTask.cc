/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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
  JSDEBUG << "Deleting task with id=" << GetId() << " and TaskNumber= " << GetMyTaskNumber();
}

void JetScapeTask::Init()
{
  JSDEBUG;
}

  /** Recursive initialization of all the subtasks of the JetScapeTask. Subtasks are also of type JetScapeTask such as Pythia Gun, Trento, Energy Loss Matter and Martini etc.
   */ 
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
  for (auto it : tasks){
    JSDEBUG << "Executing " << it->GetId();
    if (it->active_exec) it->Exec();
  }
}

void JetScapeTask::ClearTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  for (auto it : tasks)
    if (it->active_exec) it->Clear();
}

void JetScapeTask::WriteTasks(weak_ptr<JetScapeWriter> w)
{
  //VERBOSE(10);
  if (active_exec) {
    for (auto it : tasks)
      it->WriteTask(w);
  }
}

void JetScapeTask::CollectHeaders(weak_ptr<JetScapeWriter> w)
{
  //VERBOSE(10);
  if (active_exec) {
    for (auto it : tasks)
      it->CollectHeader(w);
  }
}

  void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks)
{
  tasks.push_back(m_tasks);
}

} // end namespace Jetscape
