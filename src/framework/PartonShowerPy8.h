// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

//Test PartonShower with graphh from GTL

#ifndef PARTONSHOWERPY8_H
#define PARTONSHOWERPY8_H

#include "PartonShower.h"
#include "Pythia8/Pythia.h"

// Remarks:
// Some issues in test program concerning static/dynamic casts ..
// Not the most efficient implementation! Make faster and also try to use iterators instead of vectors ...

namespace Jetscape {
  
class PartonShowerPy8 : public PartonShower
{

public:

 PartonShowerPy8() : PartonShower() {VERBOSESHOWER(8);isPy8Shower=true;isTrans=false;isClean=false;isTime=false;}
 virtual ~PartonShowerPy8();

 bool IsPythia8Shower() {return isPy8Shower;}
 bool IsTransformed() {return isTrans;}
 bool IsCleaned() {return isClean;}
 bool IsTime() {return isTime;}
 
 void FillShower(Pythia8::Event &event, int ps_id);
 
 void TransformPythiaGraph();
 void TransformPythiaGraphOld();
 void TransformAndDropPythiaGraph(double zMin=0.1, double R=100.);
 void CleanPythiaGraph();
 void CleanPythiaGraphOld();
 void AddTime(); // dummy so far, no time, just "generations" ...

 void TransformCleanAddTime() {TransformPythiaGraph();CleanPythiaGraph();AddTime();}
   
private:
  
  bool isPy8Shower;
  bool isTrans;
  bool isClean;
  bool isTime;
  
};
}
#endif
