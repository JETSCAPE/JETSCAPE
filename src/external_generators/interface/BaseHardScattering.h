#ifndef BaseHardScattering
#define BaseHardScattering

#include <HepMC/GenEvent.h>

class BaseHardScattering {

public :
  explicit BaseHardScattering(Parameter parameter_list);
  ~BaseHardScattering();

  // Pretty sure a reference is fine if we want to add particles to the event
  virtual void generatePartons(HepMC::GenEvent & event);

private:
  Parameter pset; // some copy of parameters set
//probably need to keep track of RNG state TODO
  
}

#endif
