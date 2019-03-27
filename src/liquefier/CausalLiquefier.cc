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
// -----------------------------------------
// This is a causal liquefier with the JETSCAPE framework
// -----------------------------------------

#include "CausalLiquefier.h"

namespace Jetscape {

void CausalLiquefier::smearing_kernel(
        Jetscape::real tau, Jetscape::real x, Jetscape::real y,
        Jetscape::real eta, const Droplet drop_i,
        std::array<Jetscape::real, 4> &jmu) const {
    jmu = {0., 0, 0, 0};
}

};
