// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// generic jetscape writer base clase

#ifndef JETSCAPEWRITER_H
#define JETSCAPEWRITER_H

#include <string>
#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include "JetClass.hpp"

namespace Jetscape {

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

  virtual void WriteInitFileXML() {};
  virtual void Write(weak_ptr<Parton> p) {};
  virtual void Write(weak_ptr<Jet> j) {};
  virtual void Write(weak_ptr<Vertex> v) {};
  virtual void Write(weak_ptr<PartonShower> ps) {};
  virtual void Write(string s) {};
  virtual void WriteComment(string s) {};
  virtual void WriteWhiteSpace(string s) {};
  virtual void Write(ostream *o) {};
  
  virtual void WriteEvent() {};
  // to be defined what data structure ...
  // or via passing writer to all modules
  // and handling data writing there ...
  // current approach ... 

  //virtual GetFileStream() {};
  
 private:

  string file_name_out;
  
};


} // end namespace Jetscape

#endif
