// Copyright @ Chun Shen
#ifndef TEST_HYDRO_FILE_JETSCAPE_H_
#define TEST_HYDRO_FILE_JETSCAPE_H_

#include "fluid_dynamics.h"
#include "Hydroinfo_h5.h"
#include "Hydroinfo_MUSIC.h"
#include "ParameterReader.h"

class HydroFile: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for the JETSCAPE integrated framework
 private:
    ParameterReader *paraRdr;
    int load_viscous;
    int hydro_type;
    
    double T_c;
    HydroinfoH5 *hydroinfo_h5_ptr;
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr; 

 public:
     HydroFile();
     ~HydroFile();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void get_hydro_info(real t, real x, real y, real z,
                         FluidCellInfo* fluid_cell_info_ptr);
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};

};

#endif  // TEST_HYDRO_FILE_JETSCAPE_H_
