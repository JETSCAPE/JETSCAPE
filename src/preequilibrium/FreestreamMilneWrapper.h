
#ifndef FREESTREAMMILNEWRAPPER_H
#define FREESTREAMMILNEWRAPPER_H

#include "PreequilibriumDynamics.h"
#include "FreestreamMilne.cpp" 

using namespace Jetscape;

class FreestreamMilneWrapper: public PreequilibriumDynamics {
 private:
  //int mode; //!< records running mode
  FREESTREAMMILNE *fsmilne_ptr;
  
 public:
  FreestreamMilneWrapper();
  ~FreestreamMilneWrapper();
  
  void initialize_preequilibrium(PreEquilibriumParameterFile parameter_list);  
  void evolve_preequilibrium();
};

#endif  // FREESTREAMMILNEWRAPPER_H
