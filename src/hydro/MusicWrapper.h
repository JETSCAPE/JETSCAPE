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

#ifndef MUSICWRAPPER_H
#define MUSICWRAPPER_H

#include "FluidDynamics.h"
#include "music.h"

using namespace Jetscape;

//! this is wrapper class for MUSIC so that it can be used as a external
//! library for the JETSCAPE integrated framework
class MpiMusic: public FluidDynamics {
 private:
    // int mode;            //!< records running mode
    MUSIC *music_hydro_ptr;
    int doCooperFrye;    //!< flag to run Cooper-Frye freeze-out
                         //!< for soft particles

 public:
     MpiMusic();
     ~MpiMusic();

     void InitializeHydro(Parameter parameter_list);

     void EvolveHydro();
     void GetHydroInfo(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
		std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);
     void GetHyperSurface(Jetscape::real T_cut,
                           SurfaceCellInfo* surface_list_ptr) {};
     void collect_freeze_out_surface();
};

#endif // MUSICWRAPPER_H
