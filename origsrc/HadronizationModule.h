

#ifndef HADRONIZATIONMODULE_H
#define HADRONIZATIONMODULE_H

#include "Hadronization.h"

namespace Jetscape {

template <typename Derived>
class HadronizationModule : public Hadronization
{
  
 public:

  using Hadronization::Hadronization;
  
  virtual shared_ptr<Hadronization> Clone() const override
   {
     auto ret = make_shared<Derived>(static_cast<const Derived&>(*this));
     return ret;
   }  
};


} // end namespace Jetscape


#endif
