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
//#include "JetScapeLogger.h"
#include "RealType.h"


namespace Jetscape {

LiquefierBase::LiquefierBase(std::array<Jetscape::real, 4> x_in,
                             std::array<Jetscape::real, 4> p_in) {
    xmu = x_in;
    pmu = p_in;
}

virtual void LiquefierBase::smearing_kernel(
        Jetscape::real tau, Jetscape::real x, Jetscape::real y,
        Jetscape::real eta, std::array<Jetscape::real, 4> &jmu) {
    jmu = {0, 0, 0, 0};
}

}

