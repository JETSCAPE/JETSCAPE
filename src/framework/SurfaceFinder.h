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
// This is a general basic class for a hyper-surface finder

#ifndef SURFACEFINDER_H_
#define SURFACEFINDER_H_

#include "RealType.h"
#include "SurfaceCellInfo.h"
#include "EvolutionHistory.h"

namespace Jetscape {

class SurfaceFinder {
 private:
     Jetscape::real T_cut;
     const EvolutionHistory &bulk_info;

 public:
    SurfaceFinder(const Jetscape::real T_in,
                  const EvolutionHistory &bulk_data);
    ~SurfaceFinder() {};

    //bool check_intersect(
    //        Jetscape::real tau, Jetscape::real x, Jetscape::real y,
    //        Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
    //        Jetscape::real ***cube);
    //int Find_full_hypersurface();
};

}

#endif  // SURFACEFINDER_H_
