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
  /** Default Constructor for JetScapeTask is defined. It sets the  flag "active_exec" to true. String variable "id" is not defined yet. It prints VERBOSE(9).               */
JetScapeTask::JetScapeTask()
{
  active_exec=true;
  id="";
  my_task_number_ = JetScapeTaskSupport::Instance()->RegisterTask();
  VERBOSE(9);
}

  /** Default Destructor for JetScapeTask is defined. It prints VERBOSE(9).
  */
JetScapeTask::~JetScapeTask()
{
  VERBOSE(9);
  DEBUG << "Deleting task with id=" << GetId() << " and TaskNumber= " << get_my_task_number();
}


  /** This function prints a DEBUG statement. It can be override by different modules' initialization function.
   */
void JetScapeTask::Init()
{
  DEBUG;
}

  /** Recursive initialization of all JetScapeTask-type tasks like PGun, Matter, Martini etc. It can be override by other modules to initialize their subtasks in JetScape framework. It also prints the number of subtasks in the current task.                  */
void JetScapeTask::InitTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
 
  //In short ...
  for (auto it : tasks)
    it->Init();
}

  /** It only prints VERBOSE(7). It can be override by different modules' Execution function.                                                                  
  */
void JetScapeTask::Exec()
{
  VERBOSE(7);
}

  /** Recursive Execution of all JetScapeTask-type tasks like PGun, Matter, Martini etc. It can be override by other modules to execute their subtasks in JetScape framework. It also prints the number of subtasks in the current task.
 */
void JetScapeTask::ExecuteTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  if (active_exec)
    for (auto it : tasks)
      it->Exec();
}

  /** Recursively erase  the memory and pointers allocated for the JetScapeTask-type subtasks like PGun, Matter, Martini etc. It can be override by other modules to erase their memory and the pointers allocated to the subtasks in JetScape framework. It also prints the number of subtasks in the current task.               
  */
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

  /** This function prints the final state partons for each of the tasks. To print the final state partons, it calls GetFinalPartons().
  */
void JetScapeTask::GetPartons(weak_ptr<PartonPrinter> p)
{
  //cout<<"############### Printing partons in shower " << "\n";
  for (auto it : GetTaskList())
  {
    it->GetFinalPartons(p);
  }
}

  /** This function adds different modules into the vector array of JetScapeTasks. These modules should  inherit the JetScapeTask class, and override Init(), Exec(), Clear(), and WriteTask() functions.  
  */
void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks)
{
  tasks.push_back(m_tasks);
}

} // end namespace Jetscape
