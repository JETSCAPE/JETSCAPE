// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPE_H
#define JETSCAPE_H

#include "JetScapeLogger.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeModuleBase.h"

using namespace std;

namespace Jetscape {

class JetScape : public JetScapeModuleBase
{
  
 public:
  
  JetScape();
  JetScape(string m_name) : JetScapeModuleBase (m_name) {n_events=1;VERBOSE(8);}
  JetScape(string m_name, int m_n_events) : JetScapeModuleBase (m_name) {n_events=m_n_events;VERBOSE(8);}
  virtual ~JetScape();

  void Init(); 
  void Exec();
  void Finish();

  void SetNumberOfEvents(int m_n_events) {n_events=m_n_events;}
  int GetNumberOfEvents() {return n_events;}

 private:

  void SetPointers();
  
  void Show();
  int n_events;
  
};

} // end namespace Jetscape

#endif
