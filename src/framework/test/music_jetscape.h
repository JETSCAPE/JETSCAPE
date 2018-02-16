// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef TEST_MUSIC_JETSCAPE_H_
#define TEST_MUSIC_JETSCAPE_H_

#include "FluidDynamics.h"
#include "music.h"

using namespace Jetscape;

//! this is wrapper class for MUSIC so that it can be used as a external
//! library for the JETSCAPE integrated framework
class MPI_MUSIC: public FluidDynamics {
 private:
    int mode;            //!< records running mode
    MUSIC *music_hydro_ptr;
    int doCooperFrye;    //!< flag to run Cooper-Frye freeze-out
                         //!< for soft particles

 public:
     MPI_MUSIC();
     ~MPI_MUSIC();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void get_hydro_info(real t, real x, real y, real z,
                         std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // TEST_MUSIC_JETSCAPE_H_
