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
  /** Default constructor to create the main task of the JetScape framework. It sets the total number of events to 1.
  */  
  JetScape();
  /** This is a constructor to create the main task of the JetScape framework. It sets the XML file name to "m_name", and total number of events to 1.   
   */
  JetScape(string m_name) : JetScapeModuleBase (m_name) {n_events=1;VERBOSE(8);}
  /** This is a constructor to create the main task of the JetScape framework. It sets the XML file name to "m_name", and total number of events to "m_n_events".
   */
  JetScape(string m_name, int m_n_events) : JetScapeModuleBase (m_name) {n_events=m_n_events;VERBOSE(8);}

  /** This is a destructor for a JetScape.
   */
  virtual ~JetScape();
  
  /** This function initializes the main task of the JetScape framework. It calls JetScapeTask::InitTaks() function to initialize the modules/tasks of a JetScapeTask.
  */
  void Init(); 

  /** This function execute the modules/tasks of a JetScapeTask for all the events. It also calls "GetPartons()" function to print parton shower, and  "WriteTasks()" function to store the data in the XML file.  
  */  
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
