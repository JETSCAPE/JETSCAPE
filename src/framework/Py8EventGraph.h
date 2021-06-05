// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL

#ifndef PY8EVENTGRAPH_H
#define PY8EVENTGRAPH_H

#include "JetScapeGraph.h"
#include "Pythia8/Pythia.h"

// Remarks:
// Some issues in test program concerning static/dynamic casts ..
// Not the most efficient implementation! Make faster and also try to use iterators instead of vectors ...

namespace Jetscape {

class Py8EventGraph : public JetScapeGraph<JetScapeParticleBase>
{

public:

 Py8EventGraph() : JetScapeGraph<JetScapeParticleBase>() {VERBOSESHOWER(8);isPy8Shower=true;isTrans=false;isClean=false;isTime=false;}
 virtual ~Py8EventGraph();

 bool IsPythia8Shower() {return isPy8Shower;}
 bool IsTransformed() {return isTrans;}
 bool IsCleaned() {return isClean;}
 bool IsTime() {return isTime;}

 void FillEvent(Pythia8::Event &event);

private:

  bool isPy8Shower;
  bool isTrans;
  bool isClean;
  bool isTime;

};
}
#endif
