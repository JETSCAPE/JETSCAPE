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

  /** Default constructor to create a JetScapeTask. It sets the flag "active_exec" to true and  "id" to default string value.          
   */
JetScapeTask::JetScapeTask()
{
  active_exec=true;
  id="";
  my_task_number_ = JetScapeTaskSupport::Instance()->RegisterTask();
  VERBOSE(9);
}

  /** Default destructor for a JetScapeTask.     
   */
JetScapeTask::~JetScapeTask()
{
  VERBOSE(9);
  DEBUG << "Deleting task with id=" << GetId() << " and TaskNumber= " << get_my_task_number();
}

  /**  A virtual function to define a default initialization function for a JetScapeTask. It can  be overridden by different modules/tasks.
   */
void JetScapeTask::Init()
{
  DEBUG;
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

  /** A virtual function to define a default Exec() function for a JetScapeTask. It can be overridden by different modules/tasks.
   */
void JetScapeTask::Exec()
{
  VERBOSE(7);
}

  /** Recursive Execution of all the subtasks of the JetScapeTask.
   */
void JetScapeTask::ExecuteTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  if (active_exec)
    for (auto it : tasks)
      it->Exec();
}

  /** Recursively calls Clear() function of the subtasks of a JetScapeTask.
   */
void JetScapeTask::ClearTasks()
{
  VERBOSE(7) << " : # Subtasks = "<<tasks.size();
  for (auto it : tasks)
    it->Clear();
}

  /** Recursively write the output information of different tasks/subtasks of a JetScapeTask into a file. We use "active_exec" flag  to decide whether to write the output in the file or not.
   */
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


  /** This function prints the partons shower.
   */
void JetScapeTask::GetPartons(weak_ptr<PartonPrinter> p)
{
  //cout<<"############### Printing partons in shower " << "\n";
  for (auto it : GetTaskList())
  {
    it->GetFinalPartons(p);
  }
}

  /** This function adds the module "m_tasks" into the vector of subtask of a JetScapeTask. 
   */

void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks)
{
  tasks.push_back(m_tasks);
}

} // end namespace Jetscape
