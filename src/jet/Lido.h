// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...


#ifndef LIDO_H
#define LIDO_H

#include <fstream>
#include <math.h>
#include "JetEnergyLossModule.h"
#include "workflow.h"

using namespace Jetscape;

class LidoVirtualList: public fjcore::PseudoJet::UserInfoBase
  {
  public:
    std::vector<particle> radlist_;
    LidoVirtualList(vector<particle> radlist)  :radlist_(radlist) {};
    std::vector<particle> radlist() const {return radlist_;}
  };

class Lido : public JetEnergyLossModule<Lido>
{
 private:

 std::map<int, std::vector<Process> > AllProcesses;
 double hydro_Tc;      // critical temperature
 double hydro_tStart;  // initilization time of hydro
 double Q0;

 public:

  Lido();
  virtual ~Lido();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
};



#endif
