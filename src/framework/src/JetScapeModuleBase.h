// JetScapeModuleBase class

#ifndef JETSCAPEMODULEBASE_H
#define JETSCAPEMODULEBASE_H

#include <string>
#include <memory>

#include "JetScapeTask.h"
//#incldue "JetScapeWriter.h"
#include "sigslot.h"

class JetScapeWriter;

using namespace std;

class JetScapeModuleBase : public JetScapeTask , public sigslot::has_slots<sigslot::multi_threaded_local>
{
  
 public: 

  JetScapeModuleBase();
  JetScapeModuleBase(string m_name);
  virtual ~JetScapeModuleBase();
  
  virtual void Init();
  virtual void Exec() {};
  virtual void Clear() {};

  // Think about workflow ... or in Task !!!????
  //virtual void WriteTasks(weak_ptr<JetScapeWriter> w) {};
  //virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
  //void OpenXMLFile();
  void SetXMLFileName(string m_name) { xml_file_name = m_name; }
  string GetXMLFileName() {return xml_file_name;}
  //static void SetCurrentEvent(int m_current_event) {current_event=m_current_event;}
  static int GetCurrentEvent() {return current_event;}
  static void IncrementCurrentEvent() {current_event++;}
  
 private:

  string xml_file_name;
  static int current_event; //(better static!??)
  
};

#endif
