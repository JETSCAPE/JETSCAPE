// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...

#ifndef DUMMYELOSS_H
#define DUMMYELOSS_H

#include "JetEnergyLossModule.h"

using namespace Jetscape;

class DummyEloss : public JetEnergyLossModule<DummyEloss> //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  DummyEloss();
  virtual ~DummyEloss();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {}; //funny, should not break if not not overriden !???
  
 protected:
  
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif

