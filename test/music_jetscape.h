// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef TEST_MUSIC_JETSCAPE_H_
#define TEST_MUSIC_JETSCAPE_H_

#include "music/src/util.h"
#include "music/src/grid.h"
#include "music/src/data.h"
#include "music/src/init.h"
#include "music/src/eos.h"
#include "music/src/evolve.h"
#include "../src/fluid_dynamics.h"

class MUSIC: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for the JETSCAPE integrated framework
 private:
     int mode;            // records running mode

     InitData *DATA;
     Util *util;
     
     Grid ***arena;

     EOS *eos;
     Init *init;
     Evolve *evolve;

public:
     MUSIC(Parameter parameter_list);
     ~MUSIC();

     void ReadInData3(string file);
     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void get_hydro_info(real t, real x, real y, real z,
                         FluidCellInfo* fluid_cell_info_ptr) {};
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};

};

#endif  // TEST_MUSIC_JETSCAPE_H_
