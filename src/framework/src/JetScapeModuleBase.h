// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPEMODULEBASE_H
#define JETSCAPEMODULEBASE_H

#include <string>
#include <memory>
#include <random>

#include "JetScapeTask.h"
#include "sigslot.h"

namespace Jetscape {

class JetScapeWriter;

class JetScapeModuleBase : public JetScapeTask , public sigslot::has_slots<sigslot::multi_threaded_local>
{
  
 public: 

  JetScapeModuleBase();
  JetScapeModuleBase(string m_name);
  virtual ~JetScapeModuleBase();

  //virtual shared_ptr<JetScapeModuleBase> Clone() const {return nullptr;}  
  
  virtual void Init();
  /** A virtual function to define a default Exec() function for a JetScapeModuleBase. It can be overridden by different modules/tasks.
   */
  virtual void Exec() {};
  /**  A virtual function to define a default Clear() function for a JetScapeModuleBase. It can be overridden by different modules/tasks.
   */
  virtual void Clear() {};
  /** This function sets the name of the XML file to be used to store output information for the modules/tasks of a JetScapeTask.
   */ 
  void SetXMLFileName(string m_name) { xml_file_name = m_name; }
  /** This function returns the XML file name. This file contains the output data for the modules/tasks of a JetScapeTask.
   */
  string GetXMLFileName() {return xml_file_name;}
  /** This function returns the current event number.
   */
  static int GetCurrentEvent() {return current_event;}
  /** This function increases the current event number by one.
   */
  static void IncrementCurrentEvent() {current_event++;}

  shared_ptr<std::mt19937> get_mt19937_generator();
  
 private:

  string xml_file_name;
  static int current_event; 
  shared_ptr<std::mt19937> mt19937_generator_;

  
};

} // end namespace Jetscape

#endif
