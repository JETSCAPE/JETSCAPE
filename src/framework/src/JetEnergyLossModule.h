// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// Use CRTP for cloning of derived class in base class ...

#ifndef JETENERGYLOSSMODULE_H
#define JETENERGYLOSSMODULE_H

#include "JetEnergyLoss.h"

template <typename Derived>
class JetEnergyLossModule : public JetEnergyLoss
{
  
 public:

  using JetEnergyLoss::JetEnergyLoss;
  
  virtual shared_ptr<JetEnergyLoss> Clone() const override
   {
     //compiles and seems to work (use of *this bad with shared !????)
     return make_shared<Derived>(static_cast<const Derived&>(*this));
   }
     

};

#endif
