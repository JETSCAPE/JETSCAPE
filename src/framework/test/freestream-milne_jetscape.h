// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef TEST_FREESTREAM_JETSCAPE_H_
#define TEST_FREESTREAM_JETSCAPE_H_

//#include "PreequilibriumDynamics.h"
//#include "freestream-milne.h"

using namespace Jetscape;

//! this is wrapper class for MUSIC so that it can be used as a external
//! library for the JETSCAPE integrated framework
class FREESTREAM: public PreequilibriumDynamics {
 private:
    int mode;            //!< records running mode
    //MUSIC *music_hydro_ptr;
    //int doCooperFrye;    //!< flag to run Cooper-Frye freeze-out
                         //!< for soft particles

 public:
     FREESTREAM();
     ~FREESTREAM();

     void initialize_preequilibrium(Parameter parameter_list);

     void evolve_preequilibrium();
     //void get_hydro_info(real t, real x, real y, real z, FluidCellInfo* fluid_cell_info_ptr);
     //void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // TEST_FREESTREAM_JETSCAPE_H_
