// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...

#ifndef ELOSSVALIDATION_H
#define ELOSSVALIDATION_H

#include "JetEnergyLossModule.h"

using namespace Jetscape;

class ElossValidate : public JetEnergyLossModule<ElossValidate>
{  
 public:
  
  ElossValidate();
  virtual ~ElossValidate();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);


 protected:
  
};


#endif // ELOSSVALIDATION_H

