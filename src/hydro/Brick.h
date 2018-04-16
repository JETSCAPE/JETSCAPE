/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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

     void InitializeHydro(Parameter parameter_list);

     void EvolveHydro();
     void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
			 std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void GetHypersurface(Jetscape::real T_cut, SurfaceCellInfo* surface_list_ptr) {};

     void InitTask();
     //virtual void Exec();
     
};

#endif  // BRICK_H
