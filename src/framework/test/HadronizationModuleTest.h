#ifndef HADRONIZATIONMODULESTEST_H
#define HADRONIZATIONMODULESTEST_H

#include "HadronizationModule.h"

using namespace Jetscape;

class HadronizationModuleTest : public HadronizationModule<HadronizationModuleTest>
{  
 public:
  HadronizationModuleTest();
  virtual ~HadronizationModuleTest();
  
  void Init();
  void DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);


};


#endif
