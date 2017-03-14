// generic jetscape writer base clase

#ifndef JETSCAPEWRITER_H
#define JETSCAPEWRITER_H

#include <string>
#include "JetScapeModuleBase.h"
#include "JetClass.hpp"

class JetScapeWriter : public JetScapeModuleBase
{

 public:

  JetScapeWriter() {};
  JetScapeWriter(string m_file_name_out) {file_name_out =  m_file_name_out;}
  virtual ~JetScapeWriter() {};
  
  void SetOutputFileName(string m_file_name_out) {file_name_out =  m_file_name_out;}
  string GetOutputFileName() {return file_name_out;}

  virtual bool GetStatus()=0;
  virtual void Close() {};
  virtual void Open() {};

  virtual void Write(shared_ptr<Parton> p) {};
  virtual void Write(shared_ptr<Jet> j) {};
  virtual void Write(shared_ptr<Vertex> v) {};
  virtual void Write(string s) {};
  
  virtual void WriteEvent() {};
  // to be defined what data structure ...
  // or via passing writer to all modules
  // and handling data writing there ...
  // current approach ... d
  
 private:

  string file_name_out;
  
};
#endif
