// jetscape writer ascii + gzip class

#ifndef JETSCAPEWRITERASCIIGZ_H
#define JETSCAPEWRITERASCIIGZ_H

#include "gzstream.h"

#include <fstream>
#include <string>

#include "JetScapeWriter.h"

// Maybe a templeate for Ascii and AsciiGZ !? Or in general !???

class JetScapeWriterAsciiGZ : public JetScapeWriter
{

 public:

  JetScapeWriterAsciiGZ() {};
  JetScapeWriterAsciiGZ(string m_file_name_out);
  virtual ~JetScapeWriterAsciiGZ();

  void Init();
  void Exec();
  
  bool GetStatus() {return output_file.good();}
  void Close() {output_file.close();}

  void Write(shared_ptr<Parton> p);
  void Write(string s) {output_file<<s<<endl;}
  // ...
  
  void WriteEvent();
  
 private:

  ogzstream output_file; //!< Output file
  //int m_precision; //!< Output precision
  
};
#endif
