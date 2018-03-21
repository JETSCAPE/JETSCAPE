
#ifndef TEST_FREESTREAM_JETSCAPE_H_
#define TEST_FREESTREAM_JETSCAPE_H_

#include "PreequilibriumDynamics.h"
#include "FreestreamMilne.cpp" 

using namespace Jetscape;

class FREESTREAM: public PreequilibriumDynamics {
 private:
    int mode; //!< records running mode
    FREESTREAMMILNE *fsmilne_ptr;

 public:
    FREESTREAM();
    ~FREESTREAM();

    void initialize_preequilibrium(PreEquilibriumParameterFile parameter_list);

    void evolve_preequilibrium();
};

#endif  // TEST_FREESTREAM_JETSCAPE_H_
