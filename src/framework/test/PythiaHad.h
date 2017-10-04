#ifndef PYTHIAHAD_H
#define PYTHIAHAD_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class PythiaHad : public HadronizationModule<PythiaHad>, public Pythia8::Pythia
//class PythiaHad : public HadronizationModule
{  
 public:
  
  PythiaHad(string xmlDir = "DONTUSETHIS", bool printBanner = false)
    : Pythia8::Pythia(xmlDir,printBanner), HadronizationModule()
  {
    SetId("UninitializedPythiaHad");
  }


  //PythiaHad();
  virtual ~PythiaHad();
  
  void Init();
  void DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);


};


#endif



















