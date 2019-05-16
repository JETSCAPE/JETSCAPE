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

LiquefierBase::LiquefierBase() : drop_stat(-11) {
    GetHydroCellSignalConnected = false;
}

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


void LiquefierBase::add_hydro_sources(std::vector<Parton> &pIn,
                                      std::vector<Parton> &pOut) {
    // if e_threshold > 0, use e_threshold, else, use e_threshold*T 
    // this should be put into xml later.
    auto e_threshold = 2.0;

    if (pOut.size() == 0) return;  // the process is freestreaming, ignore
    
    //cout << "debug, before ......." << pIn.size() << "  "
    //     << pOut.size() << endl;
    for (auto &iparton : pOut) {
        if (iparton.pstat() == -1) {
            // remove negative particles from parton list
            iparton.set_stat(drop_stat);
	        continue;
        }

        // for positive particles, including jet partons and recoil partons
        if (e_threshold > 0.) {
	        if (iparton.e() < e_threshold) {
                iparton.set_stat(drop_stat);
		        continue;
            }
        }
        auto tLoc = iparton.x_in().t();
        auto xLoc = iparton.x_in().x();
        auto yLoc = iparton.x_in().y();
        auto zLoc = iparton.x_in().z();
        //cout << "debug1" << endl;
        std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
        GetHydroCellSignal(tLoc, xLoc, yLoc, zLoc, check_fluid_info_ptr);
        //cout << "debug2  " << tLoc << "  " << xLoc << "  " << yLoc
        //     << "  " << zLoc << "  " << check_fluid_info_ptr << endl;
        auto tempLoc = check_fluid_info_ptr->temperature;
        //cout << "debug3" << endl;
        auto vxLoc = check_fluid_info_ptr->vx;
        auto vyLoc = check_fluid_info_ptr->vy;
        auto vzLoc = check_fluid_info_ptr->vz;
        auto beta2 = vxLoc*vxLoc + vyLoc*vyLoc + vzLoc*vzLoc;
	    auto gamma = 1.0 / sqrt(1.0 - beta2);

        // delete partons with energy smaller than 4*T
        // (in the local rest frame) from parton list
	    if (gamma*(iparton.e() - iparton.p(1)*vxLoc
                   - iparton.p(2)*vyLoc - iparton.p(3)*vzLoc) < 4.0*tempLoc) {
            iparton.set_stat(drop_stat);
		    continue;
	    }
    }
    
    const Jetscape::real hydro_source_abs_err = 1e-15;
    FourVector p_final;
    FourVector p_init;
    FourVector x_final;
    FourVector x_init;
    //cout << "debug, mid ......." << pOut.size() << endl;
    // use energy conservation to deterime the source term
    const auto weight_init = 1.0;
    for (const auto &iparton : pIn) {
        auto temp = iparton.p_in();
        p_init += temp;
        x_init = iparton.x_in();
    }
    
    auto weight_final = 0.0;
    for (const auto &iparton : pOut) {
        if (iparton.pstat() == drop_stat) continue;
        auto temp = iparton.p_in();
        p_final += temp;
        x_final = iparton.x_in(); 
        weight_final = 1.0;
    }

    if (std::abs(p_init.t() - p_final.t())/p_init.t() > hydro_source_abs_err) {
        auto droplet_t = ((x_final.t()*weight_final + x_init.t()*weight_init)
                          /(weight_final + weight_init));
        auto droplet_x = ((x_final.x()*weight_final + x_init.x()*weight_init)
                          /(weight_final + weight_init));
        auto droplet_y = ((x_final.y()*weight_final + x_init.y()*weight_init)
                          /(weight_final + weight_init));
        auto droplet_z = ((x_final.z()*weight_final + x_init.z()*weight_init)
                          /(weight_final + weight_init));

        auto droplet_tau = sqrt(droplet_t*droplet_t - droplet_z*droplet_z);
        auto droplet_eta = (0.5*log((droplet_t + droplet_z)
                                    /(droplet_t - droplet_z)));
        auto droplet_E   = p_init.t() - p_final.t(); 
        auto droplet_px  = p_init.x() - p_final.x();
        auto droplet_py  = p_init.y() - p_final.y();
        auto droplet_pz  = p_init.z() - p_final.z();

        std::array<Jetscape::real, 4> droplet_xmu = {
                            static_cast<Jetscape::real>(droplet_tau),
                            static_cast<Jetscape::real>(droplet_x),
                            static_cast<Jetscape::real>(droplet_y),
                            static_cast<Jetscape::real>(droplet_eta)};
        std::array<Jetscape::real, 4> droplet_pmu = {
                            static_cast<Jetscape::real>(droplet_E),
                            static_cast<Jetscape::real>(droplet_px),
                            static_cast<Jetscape::real>(droplet_py),
                            static_cast<Jetscape::real>(droplet_pz)};
        Droplet drop_i(droplet_xmu, droplet_pmu);
        add_a_droplet(drop_i);
    }
    //cout << "debug, after ......." << pOut.size() << endl;
}


void LiquefierBase::Clear() {
    dropletlist.clear();
}


Jetscape::real LiquefierBase::get_dropletlist_total_energy() const {
    Jetscape::real total_E = 0.0;
    for (const auto &drop_i : dropletlist) {
        total_E += drop_i.get_pmu()[0];
    }
    return(total_E);
}

};
