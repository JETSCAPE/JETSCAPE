// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...


#ifndef PARTONSHOWERGENERATOR_H
#define PARTONSHOWERGENERATOR_H

namespace Jetscape {

class JetEnergyLoss;

class PartonShowerGenerator
{
 public:

  PartonShowerGenerator() {};
  virtual ~PartonShowerGenerator() {};

  virtual void DoShower(JetEnergyLoss &j);


};

} // end namespace Jetscape

#endif
