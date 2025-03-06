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

#ifndef GLASMA_H
#define GLASMA_H

#include "PreequilibriumDynamics.h"

using namespace Jetscape;

class Glasma : public PreequilibriumDynamics {

public:
    Glasma();
    ~Glasma() {};

    void InitializePreequilibrium();

    void EvolvePreequilibrium();
    Jetscape::real GetPreequilibriumEvodtau() const {
        return dtau_;
    }

    /**
     * @return the start time of the preequilibrium stage
     * @attention In the Glasma case this is the end time, since there is no
     * preequilibrium evolution loaded to the framework. This makes the 
     * bulk_info.tau_min in the MusicWrapper.cc file think that the 
     * pre-equilibrium starts at preequilibrium_tau_max_ and then the 
     * GetHydroCellSignal function in the energy loss modules will not access 
     * the density profiles at the wrong time. 
     * Otherwise, it will store the hydro evolution starting at 
     * preequilibrium_tau_0_ and introduce a time shift in the hydro evolution.
     */
    Jetscape::real GetPreequilibriumStartTime() const {
        return preequilibrium_tau_max_;
    }

private:
    // Allows the registration of the module so that it is available to be
    // used by the Jetscape framework.
    static RegisterJetScapeModule<Glasma> reg;
    double dtau_;
};

#endif   // GLASMA_H
