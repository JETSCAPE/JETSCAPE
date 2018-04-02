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
#ifndef BRICK_H
#define BRICK_H

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
     void get_hydro_info(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
			 std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void get_hypersurface(Jetscape::real T_cut, SurfaceCellInfo* surface_list_ptr) {};

     void InitTask();
     //virtual void Exec();
     
};

#endif  // BRICK_H
