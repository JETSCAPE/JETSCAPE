// Framework test (dummy) JetEnergyLoss class (to be changed with real implemenation)

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
