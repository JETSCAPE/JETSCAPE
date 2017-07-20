// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLossModule.h"

class Matter : public JetEnergyLossModule<Matter> //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);
  
 private:

};


class Martini : public JetEnergyLossModule<Martini> //, public std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
 private:

};

#endif

