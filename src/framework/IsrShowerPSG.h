// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...


#ifndef ISRSHOWERPSG_H
#define ISRSHOWERPSG_H

#include "PartonShowerGenerator.h"
#include "PartonShower.h"
#include <GTL/edge.h>

namespace Jetscape {

class JetEnergyLoss;

class IsrShowerPSG : public PartonShowerGenerator
{
 public:

  IsrShowerPSG() : PartonShowerGenerator() {}
  virtual ~IsrShowerPSG() {};

  virtual void DoCalculateTime(JetEnergyLoss &j);
  virtual void DoExecTime(JetEnergyLoss &j);
  virtual void DoInitPerEvent(JetEnergyLoss &j);
  virtual void DoFinishPerEvent(JetEnergyLoss &j);

 private:

   //shared_ptr<PartonShower> pS = nullptr;
   //if use pointers or anything, proper copy has to be made of the PSG!!!

};

} // end namespace Jetscape

#endif
