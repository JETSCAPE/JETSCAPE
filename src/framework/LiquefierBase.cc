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
#include "LiquefierBase.h"
#include <math.h>

namespace Jetscape {

void LiquefierBase::get_source(Jetscape::real tau, Jetscape::real x,
                               Jetscape::real y, Jetscape::real eta,
                               std::array<Jetscape::real, 4> &jmu) const {
    jmu = {0.0, 0.0, 0.0, 0.0};
    for (const auto &drop_i : dropletlist) {

        const auto x_drop = drop_i.get_xmu();
        double ds2
            = tau*tau + x_drop[0]*x_drop[0]
            - 2.0*tau*x_drop[0]*cosh(eta-x_drop[3])
            - (x-x_drop[1])*(x-x_drop[1])
            - (y-x_drop[2])*(y-x_drop[2]);
        
        if( tau >= x_drop[0] && ds2 >= 0.0 ){
            std::array<Jetscape::real, 4> jmu_i = {0.0, 0.0, 0.0, 0.0};
            smearing_kernel(tau, x, y, eta, drop_i, jmu_i);
            for (int i = 0; i < 4; i++) jmu[i] += jmu_i[i];
        }
        
    }
}

};
