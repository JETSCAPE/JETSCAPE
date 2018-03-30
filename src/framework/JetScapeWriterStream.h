// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer ascii class

#ifndef JETSCAPEWRITERSTREAM_H
#define JETSCAPEWRITERSTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

namespace Jetscape {

template<class T>  
class JetScapeWriterStream : public JetScapeWriter
{

 public:

  JetScapeWriterStream<T>() {};
  JetScapeWriterStream<T>(string m_file_name_out);
  virtual ~JetScapeWriterStream<T>();

  void Init();
  void Exec();
  
  bool GetStatus() {return output_file.good();}
  void Close() {output_file.close();}

  void WriteInitFileXML();

  void Write(weak_ptr<PartonShower> ps);
  void Write(weak_ptr<Parton> p);
  void Write(weak_ptr<Vertex> v);
  void Write(weak_ptr<Hadron> h);
  void WriteHeaderToFile();
  
  void Write(string s) {output_file<<s<<endl;}
  void WriteComment(string s) {output_file<<"# "<<s<<endl;}
  void WriteWhiteSpace(string s) {output_file<<s<<" ";}
  void WriteEvent(); 
  
 protected:

  T output_file; //!< Output file
  //int m_precision; //!< Output precision
  
};

typedef JetScapeWriterStream<ofstream> JetScapeWriterAscii;
#ifdef USE_GZIP
typedef JetScapeWriterStream<ogzstream> JetScapeWriterAsciiGZ;
#endif

} // end namespace Jetscape

#endif // JETSCAPEWRITERSTREAM_H
