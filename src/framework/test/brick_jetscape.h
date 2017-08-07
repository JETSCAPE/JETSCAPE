// Copyright @ Chun Shen
#ifndef TEST_BRICK_JETSCAPE_H_
#define TEST_BRICK_JETSCAPE_H_

#include "FluidDynamics.h"

using namespace Jetscape;

class Brick: public FluidDynamics {
    // this is wrapper class for a simple brick
    // so that it can be used within the JETSCAPE framework
 private:
    double T_brick;
    double start_time;
    bool bjorken_expansion_on;

 public:
    
    Brick();
     ~Brick();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void get_hydro_info(real t, real x, real y, real z,
                         FluidCellInfo* fluid_cell_info_ptr);
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};

     void InitTask();
     //virtual void Exec();
     
};

#endif  // TEST_BRICK_JETSCAPE_H_
