/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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

  /** Controls whether to reuse a hydro event (for speedup).
      The number of times is controled by set_n_reuse_hydro
   */
  inline void set_reuse_hydro( const bool reuse_hydro ){ reuse_hydro_ = reuse_hydro;}
  /** Returns whether hydro events are reused.
   */
  inline bool get_reuse_hydro() const { return reuse_hydro_; }
  
  /** Controls number of times a hydro event gets reused.
      Reusal has to be explicitly turned on by set_reuse_hydro.
      Turn it on first to avoid a warning.
   */
  inline void set_n_reuse_hydro( const unsigned int n_reuse_hydro ){
    if ( !get_reuse_hydro() ){
      WARN << "Number of hydro reusals set, but reusal not turned on.";
      WARN << "Try jetscape->set_reuse_hydro (true);";
    }
    n_reuse_hydro_ = n_reuse_hydro;
  }
  inline unsigned int get_n_reuse_hydro() const { return n_reuse_hydro_; }

 protected:

  void SetPointers();
  
  void Show();
  int n_events;

  bool reuse_hydro_;
  unsigned int n_reuse_hydro_;
  
};

} // end namespace Jetscape

#endif
