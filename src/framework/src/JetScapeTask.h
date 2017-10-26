// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// JetScapeTask class
// Implement truly recursive (not yet done)
// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)
// (not really thogh ...)

#ifndef JETSCAPETASK_H
#define JETSCAPETASK_H

#include <list>
#include <vector>
#include <string>
#include <memory>

using std::vector;
using std::string;
using std::weak_ptr;
using std::shared_ptr;

namespace Jetscape {

// need forward declaration
class JetScapeWriter;

class PartonPrinter;
class Parton;

class JetScapeTask 
{
  
 public:

  JetScapeTask();
  virtual ~JetScapeTask();
  
  virtual void Init();
  virtual void Exec();

  /** A virtual function to define a default Finish() function for a JetScapeTask. It can  be overridden by different modules/tasks.
   */
  virtual void Finish() {};

  /** A virtual function to define a default Clear() function for a JetScapeTask. It can be overridden by different modules/tasks.
   */
  virtual void Clear() {};

  // Extensions for "recursive" handling ...
  virtual void ExecuteTasks();

  /** A virtual function to define a default ExecuteTask() function for a JetScapeTask. It can be overridden by different modules/tasks. 
   */
  virtual void ExecuteTask() {};
  /** A virtual function to define a default InitTask() function for a JetScapeTask. It can be overridden by different modules/tasks.                                 
   */
  virtual void InitTask() {};
  virtual void InitTasks();

  // really decide and think what is the best way (workflow ...)
  virtual void ClearTasks();
  /** A virtual function to define a default ClearTask() function for a JetScapeTask. It can be overridden by different modules/tasks.
   */
  virtual void ClearTask() {};
  /** A virtual function to define a default FinishTask() function for a JetScapeTask. It can be overridden by different modules/tasks.                           
   */
  virtual void FinishTask() {};
  /** A virtual function to define a default FinishTasks() function for a JetScapeTask. It can be overridden by different modules/tasks.
   */
  virtual void FinishTasks() {};
  
  //add here a write task (depending on if JetScapeWriter is initiallized and active) ...
  // Think about workflow ...
  virtual void WriteTasks(weak_ptr<JetScapeWriter> w);

  /** A virtual function to define a default WriteTask() function for a JetScapeTask. It can be overridden by different modules/tasks.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};

  virtual void Add(shared_ptr<JetScapeTask> m_tasks);
  /** This function returns the current task number. 
   */
  virtual const inline int get_my_task_number() const {return my_task_number_;} ;

  /** This function returns the vector of tasks of a JetScapeTask.
   */
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const {return tasks;}

  /** This function returns the task at  ith location in the vector of tasks of a JetScapeTask.*/
  shared_ptr<JetScapeTask> GetTaskAt(int i) {return tasks.at(i);}

  /** This function deletes the last task in the vector  of tasks of a JetScapeTask. */
  void EraseTaskLast() {tasks.erase(tasks.begin()+(tasks.size()-1));}
  //funny syntax (why last() not working here!?)

  /** This function deletes the task at ith location in the vector of tasks of a JetScapeTask.
   */
  void EraseTaskAt(int i) {tasks.erase((tasks.begin()+i));}


  /** This function resizes the length of the vector of tasks to "i". If "i" is less than the current size, it will keep the first i elements of the vector of the tasks of a JetScapeTask. 
   */
  void ResizeTaskList(int i) {tasks.resize(i);}

  /** This function removes all the tasks in the vector of tasks of a JetScapeTask and changes the size to 0.
   */
  void ClearTaskList() {tasks.clear();}

  /** This function returns the number of tasks of a JetScapeTask stored in the vector array of tasks.
   */
  int GetNumberOfTasks() {return (int) tasks.size();}

  /** This function tells whether the task is active or not.
   */
  const bool GetActive() const {return active_exec;}

  /** This functions sets the flag "active_exec" to true (active) or false(deactive).
   */
  void SetActive(bool m_active_exec) {active_exec=m_active_exec;}
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...
 
  /** This function sets the string "id" of the task of a JetScapeTask.
   */
  void SetId(string m_id) {id=m_id;}

  /** This function returns the id of the task of a JetScapeTask.
   */
  const string GetId() const {return id;}

  /** This function returns the vector of final state partons to be used as input for the recommbination/hadronization module.
   */
  vector<shared_ptr<Parton>>&  GetRecomPartons(){return fPartons;}

 private:

  // can be made sortabele to put in correct oder or via xml file ...
  vector<shared_ptr<JetScapeTask>> tasks;
  //list<shared_ptr<JetScapeTask>> tasks; // list vs vector any advantage of list?

  // final partons ready for recombination
  vector<shared_ptr<Parton>> fPartons;  

  bool active_exec;
  string id;
  // if for example a search rather position ... (or always sort with predefined order!?)

  int my_task_number_;
  
};

} // end namespace Jetscape

#endif
