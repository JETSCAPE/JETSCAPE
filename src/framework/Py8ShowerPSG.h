// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...


#ifndef PARTONSHOWERGENERATORDEFAULT_H
#define PARTONSHOWERGENERATORDEFAULT_H

#include "PartonShowerGenerator.h"
#include "PartonShower.h"
#include <GTL/edge.h>

namespace Jetscape {

class JetEnergyLoss;

class Py8ShowerPSG : public PartonShowerGenerator
{
 public:

  Py8ShowerPSG() : PartonShowerGenerator() {}
  virtual ~Py8ShowerPSG() {};

  virtual void DoShower(JetEnergyLoss &j);

 private:

  void EvolvePartonInTime(double startT, double endT, double deltaT, edge e, std::shared_ptr<PartonShower> pS, JetEnergyLoss  &j);

};

} // end namespace Jetscape

#endif
