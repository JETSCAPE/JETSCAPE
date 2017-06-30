// Copyright @ Chun Shen
#ifndef TEST_GUBSER_HYDRO_JETSCAPE_H_
#define TEST_GUBSER_HYDRO_JETSCAPE_H_

#include "JetScapeLogger.h"

#include "FluidDynamics.h"

using namespace Jetscape;

class GubserHydro: public FluidDynamics {
    // this is wrapper class for a simple brick
    // so that it can be used within the JETSCAPE framework
 private:
    double q;
    double e_0;

 public:
     GubserHydro();
     ~GubserHydro();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     double get_temperature(double e_local);
     void get_hydro_info(real t, real x, real y, real z,
                         FluidCellInfo* fluid_cell_info_ptr);
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};

};

#endif  // TEST_GUBSER_HYDRO_JETSCAPE_H_
