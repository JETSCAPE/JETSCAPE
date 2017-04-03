// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Pythia8 Test ...

#ifndef JSPYTHIA8_H
#define JSPYTHIA8_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

class JSPythia8: public HardProcess , public Pythia8::Pythia {
  
 private:
    
 public:
    
  JSPythia8();
  JSPythia8(string xmlDir = "../xmldoc", bool printBanner = false): Pythia8::Pythia(xmlDir,printBanner) , HardProcess()
    {SetId("JSPythia8");}
  ~JSPythia8();
     
  void InitTask();
  void Exec();
     
};

#endif  // JSPYTHIA8_H
