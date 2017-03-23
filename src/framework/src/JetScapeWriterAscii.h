// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer ascii class

#ifndef JETSCAPEWRITERASCII_H
#define JETSCAPEWRITERASCII_H

#include <fstream>
#include <string>

#include "JetScapeWriter.h"

class JetScapeWriterAscii : public JetScapeWriter
{

 public:

  JetScapeWriterAscii() {};
  JetScapeWriterAscii(string m_file_name_out);
  virtual ~JetScapeWriterAscii();

  void Init();
  void Exec();
  
  bool GetStatus() {return output_file.good();}
  void Close() {output_file.close();}

  void Write(shared_ptr<Parton> p);
  void Write(string s) {output_file<<s<<endl;}
  // ...
  
  void WriteEvent(); 
  std::ofstream* GetFileStream() {return &output_file;}
      
 private:

  std::ofstream output_file; //!< Output file
  //int m_precision; //!< Output precision
  
};
#endif
