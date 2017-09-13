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

using namespace std;

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
  virtual void Finish() {};
  virtual void Clear() {};

  // Extensions for "recursive" handling ...
  virtual void ExecuteTasks();
  virtual void ExecuteTask() {};
  virtual void InitTask() {};
  virtual void InitTasks();

  // really decide and think what is the best way (workflow ...)
  virtual void ClearTasks();
  virtual void ClearTask() {};
  virtual void FinishTask() {};
  virtual void FinishTasks() {};
  
  //add here a write task (depending on if JetScapeWriter is initiallized and active) ...
  // Think about workflow ...
  virtual void WriteTasks(weak_ptr<JetScapeWriter> w);
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};

  virtual void Add(shared_ptr<JetScapeTask> m_tasks);
  
  virtual const inline int get_my_task_number() const {return my_task_number_;} ;

  const vector<shared_ptr<JetScapeTask>> GetTaskList() const {return tasks;}
  shared_ptr<JetScapeTask> GetTaskAt(int i) {return tasks.at(i);}
  
  void EraseTaskLast() {tasks.erase(tasks.begin()+(tasks.size()-1));}
  //funny syntax (why last() not working here!?)
  
  void EraseTaskAt(int i) {tasks.erase((tasks.begin()+i));}
  void ResizeTaskList(int i) {tasks.resize(i);}
  void ClearTaskList() {tasks.clear();}
  
  int GetNumberOfTasks() {return (int) tasks.size();}
  
  const bool GetActive() const {return active_exec;}
  void SetActive(bool m_active_exec) {active_exec=m_active_exec;}
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...
 
  void SetId(string m_id) {id=m_id;}
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
