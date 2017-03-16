// jetscape writer ascii + gzip class

#ifndef JETSCAPEWRITERHEPMC_H
#define JETSCAPEWRITERHEPMC_H

#include <fstream>
#include <string>

#include "JetScapeWriter.h"

#include "HepMC/GenEvent.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/WriterAscii.h"
#include "HepMC/Print.h"

using namespace HepMC;

class JetScapeWriterHepMC : public JetScapeWriter , public WriterAscii
{

 public:

 JetScapeWriterHepMC() : WriterAscii("") {};
 JetScapeWriterHepMC(string m_file_name_out) : JetScapeWriter(m_file_name_out), WriterAscii(m_file_name_out) {};
  virtual ~JetScapeWriterHepMC();

  void Init();
  void Exec();
  
  bool GetStatus() {return failed();}
  void Close() {close();}

  // overload write functions ...
  void WriteEvent(); 
  
 private:

  //int m_precision; //!< Output precision
  
};
#endif
