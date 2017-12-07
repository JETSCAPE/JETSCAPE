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
