// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef SRC_MUSIC_H_
#define SRC_MUSIC_H_

#include "music/src/util.h"
#include "music/src/grid.h"
#include "music/src/data.h"
#include "music/src/init.h"
#include "music/src/eos.h"
#include "music/src/evolve.h"
#include "../src/fluid_dynamics.h"

class MUSIC: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for integrated framework, such as JETSCAPE
 private:
     int mode;            // records running mode

     InitData *DATA;
     Util *util;
     
     Grid ***arena;

     EOS *eos;
     Init *init;
     Evolve *evolve;

 public:
     MUSIC();
     ~MUSIC();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void ReadInData3(string file);
};

#endif  // SRC_MUSIC_H_
