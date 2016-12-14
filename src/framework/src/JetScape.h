// JetScape class

#ifndef JETSCAPE_H
#define JETSCAPE_H

#include "JetScapeLogger.h"
#include "JetScapeModuleBase.h"

//#include<iostream>
//#include<sstream>

using namespace std;

class JetScape : public JetScapeModuleBase
{
  
 public:
  
  JetScape();
  JetScape(string m_name) : JetScapeModuleBase (m_name) {n_events=1;VERBOSE(8);}
  JetScape(string m_name, int m_n_events) : JetScapeModuleBase (m_name) {n_events=m_n_events;VERBOSE(8);}
  virtual ~JetScape();

  void Init(); 
  void Exec();

  void SetNumberOfEvents(int m_n_events) {n_events=m_n_events;}
  int GetNumberOfEvents() {return n_events;}
  
  // dummy here, but might be helpful to implement in most classes
  // in which output is required and being conistent with C++ cout
  // vectors, fulidcells ...
  //friend ostream &operator<<(ostream &out, const JetScape &js) {out<<"Test"; return out;} 

 private:

  void SetPointers();
  
  void Show();
  int n_events;
  
};

#endif
