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
#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>

#include<iostream>
#include <sstream>

using namespace std;

namespace Jetscape {

void IsrShowerPSG::GetFinalEdgesForTime(shared_ptr<PartonShower> pS, double t, vector<edge> &vE)
{
  graph::edge_iterator eIt, eEnd;
  for (eIt = pS->edges_begin(), eEnd = pS->edges_end(); eIt != eEnd; ++eIt)
  {
    if (eIt->target().outdeg()<1)
    {
      node nEnd = eIt->target();
      auto vEnd = pS->GetVertex(nEnd);
      double tEnd = vEnd->x_in().t();
      //DEBUG
      //cout<<tEnd<<endl;
      //cout<<nE<<endl;
      if (t>tEnd)
        vE.push_back(*eIt);
    }
  }
}

// not the most efficient way via GetFinalEdgesForTime ...
void IsrShowerPSG::GetFinalPartonsForTime(shared_ptr<PartonShower> pS, double t, vector<std::shared_ptr<Parton>> &vP)
{
  vector<edge> vecE;
  GetFinalEdgesForTime(pS,t,vecE);

  for (auto e : vecE) vP.push_back(pS->GetParton(e));

  vecE.clear();
}

void IsrShowerPSG::DoCalculateTime(JetEnergyLoss &j)
{
  VERBOSE(3);
}

//REMARK: Not the most elegant way to reuse the standard DoExecTime() in JetEnergyLoss ...
//        but seems to work. Think about how to make it more efficient and avoid making things public ... !!!!
void IsrShowerPSG::DoExecTime(JetEnergyLoss &j)
{
  double currentTime = j.GetModuleCurrentTime()+j.GetModuleDeltaT();

  VERBOSE(2)<<" t = "<<currentTime;

  auto pS=j.GetShower();

  vector<edge> vecE;
  GetFinalEdgesForTime(pS,currentTime,vecE);

  if (vecE.size()>0)
  {
    for (auto e : vecE)
    {
      //cout<<e<<" ";

      j.pIn.push_back(*pS->GetParton(e));
      j.vStartVec.push_back(e.target());
    }

    //cout<<endl;
  }

  j.DoExecTime();

  j.pIn.clear(); j.vStartVec.clear();

  vecE.clear();
}

void IsrShowerPSG::DoInitPerEvent(JetEnergyLoss &j)
{
  VERBOSE(2);

  j.foundchangedorig = true;

  //debug:
  /*
  auto pS=j.GetShower();
  pS->SaveAsGV("isr_"+std::to_string(j.GetMyTaskNumber())+".gv");
  */
  //cout<<"IsrShowerPSG::DoInitPerEvent(JetEnergyLoss &j)"<<endl;
  //pS->PrintNodes(false);
  //pS->PrintEdges(false);
  //cout<<this<<endl;

}

// for debug ...
/*
void IsrShowerPSG::DoFinishPerEvent(JetEnergyLoss &j)
{
  VERBOSE(2);

  auto pS=j.GetShower();
  pS->SaveAsGV("isr_fsr_"+std::to_string(j.GetMyTaskNumber())+".gv");

  //pS->PrintEdges(false);
  //cout<<&j<<endl;
}
*/

} // end namespace Jetscape
