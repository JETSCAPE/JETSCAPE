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

  /** It does nothing. This function can be override by the different modules Finish() function.               
  */
  virtual void Finish() {};

  /** It does nothing. This function can be override by the different modules Clear() function.                                  
  */
  virtual void Clear() {};

  // Extensions for "recursive" handling ...
  virtual void ExecuteTasks();

  /** It does nothing. This function can be override by the different modules ExecuteTask() function.                            
  */
  virtual void ExecuteTask() {};
  /** It does nothing. This function can be override by the different modules InitTask() function.
  */
  virtual void InitTask() {};
  virtual void InitTasks();

  // really decide and think what is the best way (workflow ...)
  virtual void ClearTasks();
  /** It does nothing. This function can be override by the different modules ClearTask() function.
  */
  virtual void ClearTask() {};
  /** It does nothing. This function can be override by the different modules FinishTask() function.                             
  */
  virtual void FinishTask() {};
  /** It does nothing. This function can be override by the different modules FinishTasks() function.                             
  */
  virtual void FinishTasks() {};
  
  //add here a write task (depending on if JetScapeWriter is initiallized and active) ...
  // Think about workflow ...
  virtual void WriteTasks(weak_ptr<JetScapeWriter> w);
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};

  // Printer method that prints the partons of the shower
  // Is it only for EnergyLoss?
  virtual void GetPartons(weak_ptr<PartonPrinter> p);
  /** This function does nothing. This function can be override by the different modules GetFinalPartons() function.
   */
  virtual void GetFinalPartons(weak_ptr<PartonPrinter> p){};

  virtual void Add(shared_ptr<JetScapeTask> m_tasks);
  
  virtual const inline int get_my_task_number() const {return my_task_number_;} ;
  /** This function returns the vector array of tasks.
   */
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const {return tasks;}
  /** This function returns the task at index i from vector array of tasks.
   */
  shared_ptr<JetScapeTask> GetTaskAt(int i) {return tasks.at(i);}
  /** This function deletes the last task in the vector array of the tasks. 
   */
  void EraseTaskLast() {tasks.erase(tasks.begin()+(tasks.size()-1));}
  //funny syntax (why last() not working here!?)

  /** This function deletes the task at index i in the vector array of the tasks. 
   */  
  void EraseTaskAt(int i) {tasks.erase((tasks.begin()+i));}

  /** This function resizes the length of the vector array of tasks to i. If i is less than the current size, it will keep the first i elements of the vector array of the tasks.
  */
  void ResizeTaskList(int i) {tasks.resize(i);}
  /** This function removes all the tasks in the vector array of tasks and changes the size to 0. 
  */
  void ClearTaskList() {tasks.clear();}
  /** This function returns the number of tasks stored in the vector array of tasks.  */
  int GetNumberOfTasks() {return (int) tasks.size();}

  /** This function returns the status of the "active_exec" flag. 
  */  
  const bool GetActive() const {return active_exec;}

  /** This functions sets the flag "active_exec" to true or false value. This flag controls both the execution of the tasks and writing tasks results into a file.
  */
  void SetActive(bool m_active_exec) {active_exec=m_active_exec;}
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...
 
  /** This function sets the ID. Not yet defined in the framework.
  */
  void SetId(string m_id) {id=m_id;}
  /** This function returns the ID. Not yet defined.             
   */
  const string GetId() const {return id;}

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
