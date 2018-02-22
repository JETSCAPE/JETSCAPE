// Copyright @ Chun Shen
#ifndef TEST_HYDRO_FILE_JETSCAPE_H_
#define TEST_HYDRO_FILE_JETSCAPE_H_

#include "FluidDynamics.h"
#include "Hydroinfo_MUSIC.h"

#include <string>

#ifdef USE_HDF5
#include "Hydroinfo_h5.h"
#endif

using namespace Jetscape;

class HydroFile: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for the JETSCAPE integrated framework
 private:
    // ParameterReader *paraRdr;
    bool load_viscous;
    int hydro_type;

    double T_c;
#ifdef USE_HDF5
    HydroinfoH5 *hydroinfo_h5_ptr;
#endif
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;

 public:
     HydroFile();
     ~HydroFile();

     //! clean up hydro event
     void clean_hydro_event();

     //! This function initials hydro parameters and read in a hydro event
     void initialize_hydro(Parameter parameter_list);


     //! This function load a VISHNew hydro event
     void read_in_hydro_event(string VISH_filename, int buffer_size,
                              int load_viscous);

     //! This function load a MUSIC hydro event
     void read_in_hydro_event(string input_file, string hydro_ideal_file,
                              int nskip_tau);

     //! This function is a dummy function
     void evolve_hydro();

     //! This function provide fluid cell information at a given
     //! space-time point
     void get_hydro_info(real t, real x, real y, real z,
                         std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // TEST_HYDRO_FILE_JETSCAPE_H_
