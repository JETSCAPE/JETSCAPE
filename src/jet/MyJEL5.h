#ifndef MyJEL5_H
#define MyJEL5_H

#include "JetEnergyLossModule.h"
#include "ElasticCollision.h"

using namespace Jetscape;

class MyJEL5 : public JetEnergyLossModule<MyJEL5>
{  
 public:
  
  MyJEL5();
  virtual ~MyJEL5();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2,
                    vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  ElasticCollision elasticCollision;

 private:
  // Allows the registration of the module so that it is available
  // to be used by the Jetscape framework.
  static RegisterJetScapeModule<MyJEL5> reg;
  
};

#endif // MyJEL5
