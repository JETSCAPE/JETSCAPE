/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// Copyright @ Chun Shen
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

 public:
     HydroFromFile();
     ~HydroFromFile();

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
     void get_hydro_info(Jetscape::real t, Jetscape::real x,
                         Jetscape::real y, Jetscape::real z,
                         std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void set_hydro_event_idx(int idx_in) {hydro_event_idx_ = idx_in;};
     int get_hydro_event_idx() {return(hydro_event_idx_);};
     void get_hypersurface(Jetscape::real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // HYDROFROMFILE_H
