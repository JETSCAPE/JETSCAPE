#ifndef PYTHIAHAD_H
#define PYTHIAHAD_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class PythiaHad : public HadronizationModule<PythiaHad>
{  
 public:

  PythiaHad();
  virtual ~PythiaHad();
  
  void Init();
  void DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

};


#endif



















