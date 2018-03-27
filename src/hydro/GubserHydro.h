// Copyright @ Chun Shen
#ifndef GUBSERHYDRO_H
#define GUBSERHYDRO_H

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

     void get_hydro_info(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
			 std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void get_hypersurface(Jetscape::real T_cut, SurfaceCellInfo* surface_list_ptr) {};

};

#endif  // GUBSERHYDRO_H
