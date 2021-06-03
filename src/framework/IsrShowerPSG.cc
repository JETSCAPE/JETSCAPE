// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...


#include "IsrShowerPSG.h"
#include "PartonShower.h"
#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include <GTL/bfs.h>

#include<iostream>
#include <sstream>

using namespace std;

namespace Jetscape {

void IsrShowerPSG::DoCalculateTime(JetEnergyLoss &j)
{
  VERBOSE(3);
}

void IsrShowerPSG::DoExecTime(JetEnergyLoss &j)
{
  VERBOSE(3);
}

void IsrShowerPSG::DoInitPerEvent(JetEnergyLoss &j)
{
  VERBOSE(2);
  auto pS=j.GetShower();

  pS->SaveAsGV("isr.gv");

  cout<<"IsrShowerPSG::DoInitPerEvent(JetEnergyLoss &j)"<<endl;
  pS->PrintNodes(false);
  //cout<<this<<endl;
}

void IsrShowerPSG::DoFinishPerEvent(JetEnergyLoss &j)
{
  VERBOSE(2);

  auto pS=j.GetShower();
  pS->PrintEdges(false);
  //cout<<&j<<endl;

  pS->SaveAsGV("isr_fsr.gv");
}

} // end namespace Jetscape
