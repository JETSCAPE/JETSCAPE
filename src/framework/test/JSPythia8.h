//Pythia8 Test ...

#ifndef JSPYTHIA8_H
#define JSPYTHIA8_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia.h"

class JSPythia8: public HardProcess , public Pythia8::Pythia {
    // this is wrapper class for a simple brick
    // so that it can be used within the JETSCAPE framework
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
