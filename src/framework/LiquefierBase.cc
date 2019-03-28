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

namespace Jetscape {

void LiquefierBase::get_source(Jetscape::real tau, Jetscape::real x,
                               Jetscape::real y, Jetscape::real eta,
                               std::array<Jetscape::real, 4> &jmu) const {
    jmu = {0.0, 0.0, 0.0, 0.0};
    for (const auto &drop_i : dropletlist) {
          std::array<Jetscape::real, 4> jmu_i = {0.0, 0.0, 0.0, 0.0};
          smearing_kernel(tau, x, y, eta, drop_i, jmu_i);
          for (int i = 0; i < 4; i++) jmu[i] += jmu_i[i];
    }
}

};
