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

using std::abs;
using std::uniform_real_distribution;



namespace Jetscape {

template <typename Derived>
class JetEnergyLossModule : public JetEnergyLoss
{
  
 public:

  using JetEnergyLoss::JetEnergyLoss;
  
  virtual shared_ptr<JetEnergyLoss> Clone() const override
   {
     DEBUG << "Cloning task with id=" << GetId() << " and TaskNumber= " << get_my_task_number();
     // DEBUG/TODO: KK: Joern's plan was to not have to call Init again, but I'm not sure that can work/is desirable.
     auto ret=make_shared<Derived>(static_cast<const Derived&>(*this));
     //ret->Init();
     return ret;
     //compiles and seems to work (use of *this bad with shared !????)
     // return make_shared<Derived>(static_cast<const Derived&>(*this));
   }
     
  // override deactivation
  void SetActive(bool m_active_exec) {
    throw std::runtime_error("SetActive not supported for energy loss modules. Please remove the module from the manager.");
  };

};

} // end namespace Jetscape

#endif
