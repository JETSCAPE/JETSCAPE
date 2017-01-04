// JetScapeTask class
// Implement truly recursive (not yet done)
// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)
// (not really thogh ...)

#ifndef JETSCAPETASK_H
#define JETSCAPETASK_H

#include <list>
#include <vector>
#include <string>

using namespace std;

class JetScapeTask 
{
  
 public:

  JetScapeTask();
  virtual ~JetScapeTask();
  
  virtual void Init();
  virtual void Exec();
  virtual void Finish() {};
  virtual void Clear() {};

  // Extensions for clearly recursive handling ...
  // Hmm, not really ad intended ...
  virtual void ExecuteTasks();
  virtual void ExecuteTask() {}; // To be seen if needed ...
  virtual void InitTasks();
  //virtual void ClearTasks() {};
  
  // Think harder and maybe in general/decide on shared vs. unique vs. raw pointer usage ...
  // Current: Just use make_shared always (propably not the most efficient solution ...)
  // Also usage here ditactes in main ...
  virtual void Add(shared_ptr<JetScapeTask> m_tasks);
  //virtual void Add(JetScapeTask *m_tasks) {tasks.push_back(m_tasks);}

  // check const wrt to deep copy constructor JetEnergyLoss
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const {return tasks;}
  //vector<shared_ptr<JetScapeTask>> GetTaskList() {return tasks;}
  shared_ptr<JetScapeTask> GetTaskAt(int i) {return tasks.at(i);}
  
  void EraseTaskLast() {tasks.erase(tasks.begin()+(tasks.size()-1));} //funny syntax (why last() not working here!?)
  void EraseTaskAt(int i) {tasks.erase((tasks.begin()+i));}
  void ResizeTaskList(int i) {tasks.resize(i);}
  void ClearTaskList() {tasks.clear();}
  
  int GetNumberOfTasks() {return (int) tasks.size();}
  const bool GetActive() const {return active_exec;}
  void SetActive(bool m_active_exec) {active_exec=m_active_exec;}
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...
  // (check if always initilazed true by default)
  
  void SetId(string m_id) {id=m_id;}
  const string GetId() const {return id;}
  
 private:

  // can be made sortabele to put in correct oder or via xml file ...
  vector<shared_ptr<JetScapeTask>> tasks;
  //list<JetScapeTask*> tasks;
  bool active_exec;
  string id; // if for example a search rather position ... (or always sort with predefined order!?)
  
};

#endif
