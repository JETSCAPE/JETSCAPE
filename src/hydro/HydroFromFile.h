/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef HYDROFROMFILE_H
#define HYDROFROMFILE_H

#include "FluidDynamics.h"
#include "Hydroinfo_MUSIC.h"

#include <string>

#ifdef USE_HDF5
#include "Hydroinfo_h5.h"
#endif

using namespace Jetscape;

class HydroFromFile: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for the JETSCAPE integrated framework
 private:
    tinyxml2::XMLElement *para_;

    int flag_read_in_multiple_hydro_;
    int hydro_event_idx_;

    int load_viscous_;
    int hydro_type_;

    int nskip_tau_;
    double T_c_;

#ifdef USE_HDF5
    HydroinfoH5 *hydroinfo_h5_ptr;
#endif
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
  
    // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<HydroFromFile> reg;

 public:
     HydroFromFile();
     ~HydroFromFile();

     //! clean up hydro event
     void clean_hydro_event();

     //! This function initials hydro parameters and read in a hydro event
     void InitializeHydro(Parameter parameter_list);


     //! This function load a VISHNew hydro event
     void read_in_hydro_event(string VISH_filename, int buffer_size,
                              int load_viscous);

     //! This function load a MUSIC hydro event
     void read_in_hydro_event(string input_file, string hydro_ideal_file,
                              int nskip_tau);

     //! This function is a dummy function
     void EvolveHydro();

     //! This function provide fluid cell information at a given
     //! space-time point
     void GetHydroInfo(Jetscape::real t, Jetscape::real x,
                         Jetscape::real y, Jetscape::real z,
                         std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void set_hydro_event_idx(int idx_in) {hydro_event_idx_ = idx_in;};
     int get_hydro_event_idx() {return(hydro_event_idx_);};
     void GetHyperSurface(Jetscape::real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // HYDROFROMFILE_H
