// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL
#include "Py8EventGraph.h"
#include <GTL/dfs.h>
#include <GTL/node.h>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>

#include "MakeUniqueHelper.h"

using namespace std;

namespace Jetscape {

// quick and dirty fix ...
//default_random_engine generator;
//uniform_real_distribution<double> distribution(0.2,1.);

Py8EventGraph::~Py8EventGraph()
{
  VERBOSESHOWER(8);//pFinal.clear();//pVec.clear();vVec.clear();
}


void Py8EventGraph::FillEvent(Pythia8::Event &event)
{
  //REMARK, not really sure, in any case, EventGraph
  // more meant for self created shower and hadronzition
  // to be explored ...
  // Something about mother and daughterlists not
  // really obvious correct in graph representation ..

  JSINFO<<" Fill Pythia Event Graph ...";

  PrintNodes();
  PrintEdges();

  vector<node> pyNodes;

  for (unsigned int ipart=0; ipart < event.size(); ++ipart)
  {
    pyNodes.push_back(new_vertex(make_shared<Vertex>(0,0,0,event[ipart].status())));
  }

  for (unsigned int ipart=1; ipart < event.size(); ++ipart)
  {
    Pythia8::Particle p = event[ipart];

    vector<int> dl=event[ipart].daughterList();
    vector<int> ml=event[ipart].motherList();

    for (int i=0;i<ml.size();i++)
    {
      if (p.isFinal() || p.isHadron())
      {
        new_particle(pyNodes[ml[i]],pyNodes[ipart],make_shared<Hadron>(event[ipart].index(), p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));
      }
      else
      {
        if (abs(p.id())==2101 || abs(p.id())==2203 || abs(p.id())==2103)
        new_particle(pyNodes[ml[i]],pyNodes[ipart],make_shared<Hadron>(event[ipart].index(), p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));
        else
        new_particle(pyNodes[ml[i]],pyNodes[ipart],make_shared<Parton>(event[ipart].index(), p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));
      }
    }

  }

  for (unsigned int ipart=0; ipart < event.size(); ++ipart)
  {
    Pythia8::Particle p = event[ipart];
    if (p.isFinal() && p.motherList().size()>1)
    {
      node nEnd=new_vertex(make_shared<Vertex>(0,0,0,event[ipart].status()));
      new_particle(pyNodes[ipart],nEnd,make_shared<Hadron>(event[ipart].index(), p.id(),p.status(),p.pT(),p.y(),p.phi(),p.e()));
    }
  }

  hide_node(pyNodes[0]);

}

}
