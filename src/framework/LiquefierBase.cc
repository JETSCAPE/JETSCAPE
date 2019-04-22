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

void LiquefierBase::add_hydro_sources(std::vector<Parton> &pIn, std::vector<Parton> &pOut) {

    double e_threshold = 2.0; // this should be put into xml later. if e_threshold > 0, use e_threshold, else, use e_threshold*T 

    if (pOut.size() == 0) return;
    //if (weak_ptr_is_uninitialized(liquefier_ptr)) return;

    //cout << "debug, before ......." << pIn.size() << "  " << pOut.size() << endl;

    const Jetscape::real hydro_source_abs_err = 1e-15;
    FourVector p_final;
    FourVector p_init;
    FourVector x_final;
    FourVector x_init;
    double weight_final, weight_init;


    int i = 0;
    while(i < pOut.size()){
        if(pOut[i].pstat() == -1){ // remove negative particles from parton list
	    pOut.erase(pOut.begin() + i);
	    continue;
	} else { // for positive particles, including jet partons and recoil partons
 
	    //double tLoc = pOut[i].x_in().t();
	    //double xLoc = pOut[i].x_in().x();
	    //double yLoc = pOut[i].x_in().y();
	    //double zLoc = pOut[i].x_in().z();
	    //cout << "debug1" << endl;
	    //std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
            //GetHydroCellSignal(tLoc, xLoc, yLoc, zLoc, check_fluid_info_ptr);
	    //cout << "debug2  " << tLoc << "  " << xLoc << "  " << yLoc << "  " << zLoc << "  " << check_fluid_info_ptr << endl;
            //double tempLoc = check_fluid_info_ptr->temperature;
	    //cout << "debug3" << endl;
            //double vxLoc = check_fluid_info_ptr->vx;
            //double vyLoc = check_fluid_info_ptr->vy;
            //double vzLoc = check_fluid_info_ptr->vz;
            //double beta2 = vxLoc*vxLoc + vyLoc*vyLoc + vzLoc*vzLoc;
	    //double gamma = 1.0 / sqrt(1.0 - beta2);

	
            //// delete partons with energy smaller than 4*T (in the local rest frame) from parton list
	    //if( gamma * (pOut[i].e() - pOut[i].p(1)*vxLoc - pOut[i].p(2)*vyLoc - pOut[i].p(3)*vzLoc) < 4.0 * tempLoc ) {
	    ////if( gamma * (pOut[i].e() - pOut[i].p_in().x()*vxLoc - pOut[i].p_in().y()*vyLoc - pOut[i].p_in().z()*vzLoc) < 4.0 * tempLoc ) {
	    if( pOut[i].e() < e_threshold ) {
		//cout << "check before remove: " << i << "  " << pOut[i].e() << "  " << pOut.size() << endl;
                pOut.erase(pOut.begin() + i);
		//cout << "check after remove: " << i << "  " << pOut.size() << endl;
		continue;
	    }
	}

	i++;
    }

    //cout << "debug, mid ......." << pOut.size() << endl;

    // use energy conservation to deterime the source term
    for (const auto &iparton : pIn) {
        auto temp = iparton.p_in();
        p_init += temp;
        x_init = iparton.x_in();
    }
    for (const auto &iparton : pOut) {
        auto temp = iparton.p_in();
	p_final += temp;
        x_final = iparton.x_in(); 
    }
    if (std::abs(p_init.t() - p_final.t())/p_init.t() > hydro_source_abs_err) { 

    	weight_init = 1.0;
	if(pOut.size() > 0) weight_final = 1.0;
	else weight_final = 0.0;
        Jetscape::real droplet_t   = (x_final.t() * weight_final + x_init.t() * weight_init) / (weight_final + weight_init);
        Jetscape::real droplet_x   = (x_final.x() * weight_final + x_init.x() * weight_init) / (weight_final + weight_init);
        Jetscape::real droplet_y   = (x_final.y() * weight_final + x_init.y() * weight_init) / (weight_final + weight_init);
        Jetscape::real droplet_z   = (x_final.z() * weight_final + x_init.z() * weight_init) / (weight_final + weight_init);
        Jetscape::real droplet_tau = (
                    sqrt(droplet_t*droplet_t - droplet_z*droplet_z));
        Jetscape::real droplet_eta = (
                    0.5*log((droplet_t + droplet_z)
                            /(droplet_t - droplet_z)));
        Jetscape::real droplet_E   = p_init.t() - p_final.t(); 
        Jetscape::real droplet_px  = p_init.x() - p_final.x();
        Jetscape::real droplet_py  = p_init.y() - p_final.y();
        Jetscape::real droplet_pz  = p_init.z() - p_final.z();

        std::array<Jetscape::real, 4> droplet_xmu = {
            droplet_tau, droplet_x, droplet_y, droplet_eta};
        std::array<Jetscape::real, 4> droplet_pmu = {
            droplet_E, droplet_px, droplet_py, droplet_pz};
        Droplet drop_i(droplet_xmu, droplet_pmu);
        //liquefier_ptr.lock()->add_a_droplet(drop_i);
        add_a_droplet(drop_i);
    }

    //cout << "debug, after ......." << pOut.size() << endl;
}


};
