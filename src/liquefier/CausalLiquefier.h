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

#ifndef CAUSALLIQUEFIER_H
#define CAUSALLIQUEFIER_H

#include "LiquefierBase.h"
#include "RealType.h"
#include <array>

#include <cmath>
#include <gsl/gsl_sf_bessel.h>

namespace Jetscape {

class CausalLiquefier: public Jetscape::LiquefierBase {
 private:

    //parameters (to be moved to xml)---------------------------
    double tau_delay;// in [fm]
    double dtau; //in [fm]

    double time_relax;// in [fm]
    double d_diff;// in [fm]
    
    double width_delta;// in [fm]
    //---------------------------
    double c_diff;
    double gamma_relax;
    
 public:
    CausalLiquefier();
    ~CausalLiquefier() {};

    void smearing_kernel(Jetscape::real tau, Jetscape::real x,
                         Jetscape::real y, Jetscape::real eta,
                         const Droplet drop_i,
                         std::array<Jetscape::real, 4> &jmu) const;

    double causal_diffusion_kernel(double t, double r) const;
    double causal_diffusion_smooth(double t, double r) const;
    double causal_diffusion_delta(double t, double r) const;

    
    double get_t(double tau, double eta) const;
    double get_z(double tau, double eta) const;
    
    double get_ptau(double px, double pz, double eta) const;
    double get_peta(double px, double pz, double eta) const;

    
};

};

#endif  // CAUSALLIQUEFIER_H
