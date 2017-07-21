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

#include "JetScapeTask.h"
#include "sigslot.h"

using namespace std;

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
  virtual void Exec() {};
  virtual void Clear() {};
  
  void SetXMLFileName(string m_name) { xml_file_name = m_name; }
  string GetXMLFileName() {return xml_file_name;}
  
  static int GetCurrentEvent() {return current_event;}
  static void IncrementCurrentEvent() {current_event++;}
  
 private:

  string xml_file_name;
  static int current_event; 
  
};

} // end namespace Jetscape

#endif
