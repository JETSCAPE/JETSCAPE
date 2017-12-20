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

namespace Jetscape {

class JetScape : public JetScapeModuleBase
{
  
 public:
  
  JetScape();
  /** This is a constructor to create the main task of the JetScape framework. It sets the XML file name to "m_name", and total number of events to 1.   
   */
  JetScape(string m_name) : JetScapeModuleBase (m_name) {n_events=1;VERBOSE(8);}
  /** This is a constructor to create the main task of the JetScape framework. It sets the XML file name to "m_name", and total number of events to "m_n_events".
   */
  JetScape(string m_name, int m_n_events) : JetScapeModuleBase (m_name) {n_events=m_n_events;VERBOSE(8);}
  virtual ~JetScape();

  void Init(); 
  void Exec();
  void Finish();

  /** This function sets the total number of events to "m_n_events".
   */
  void SetNumberOfEvents(int m_n_events) {n_events=m_n_events;}
  /** This function returns the total number of events.
   */
  int GetNumberOfEvents() {return n_events;}

 private:

  void SetPointers();
  
  void Show();
  int n_events;
  
};

} // end namespace Jetscape

#endif
