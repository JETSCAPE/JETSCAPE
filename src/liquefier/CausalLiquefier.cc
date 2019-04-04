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
#include "JetScapeLogger.h"


namespace Jetscape {
    
    
CausalLiquefier::CausalLiquefier(){
    
    //parameters (to be moved to xml)---------------------------
    tau_delay = 0.5;// in [fm]
    dtau = 0.2; //in [fm]
    
    time_relax = 0.1;// in [fm]
    d_diff = 0.08;// in [fm]
    
    width_delta = 1.0;// in [fm]
    //---------------------------
    
    c_diff = sqrt(d_diff/time_relax);
    gamma_relax = 0.5/time_relax;
}

//CausalLiquefier::~CausalLiquefier(){};
    
    
void CausalLiquefier::smearing_kernel(
        Jetscape::real tau, Jetscape::real x, Jetscape::real y,
        Jetscape::real eta, const Droplet drop_i,
        std::array<Jetscape::real, 4> &jmu) const {

    jmu = {0., 0, 0, 0};
    
    auto x_drop = drop_i.get_xmu();
    const auto p_drop = drop_i.get_pmu();
    
    double tau_drop = x_drop[0];
    double eta_drop = x_drop[3];
    
    double t = get_t(tau, eta);
    x_drop[0] = get_t(tau_drop, eta_drop);
    
    double z = get_z(tau, eta);
    x_drop[3] = get_z(tau_drop, eta_drop);

    if( tau - 0.5*dtau <= tau_drop + tau_delay &&
       tau + 0.5*dtau > tau_drop + tau_delay ){
        
        double t_delay = get_t(tau_delay, eta);
        double delta_r2 = (x-x_drop[1])*(x-x_drop[1]) + (y-x_drop[2])*(y-x_drop[2]) +(z-x_drop[3])*(z-x_drop[3]);
        
        double kernel
        = causal_diffusion_kernel(t_delay, sqrt(delta_r2));
        
        jmu[0] = kernel*get_ptau(p_drop[0], p_drop[3], eta);
        jmu[1] = kernel*p_drop[1];
        jmu[2] = kernel*p_drop[2];
        jmu[3] = kernel*get_peta(p_drop[0], p_drop[3], eta);
        
    }
    
}
    
    
double CausalLiquefier::causal_diffusion_kernel(
        double t, double r) const {
    
    double smooth = causal_diffusion_smooth(t, r);
    double delta = causal_diffusion_delta(t, r);
    
    return smooth+delta;
}
    
double CausalLiquefier::causal_diffusion_smooth(double t, double r) const {
    

    if( r < c_diff*t ){
        
        
        double u = sqrt( c_diff*c_diff*t*t - r*r );
        double x = gamma_relax*u/c_diff; // unitless

        double i1 = gsl_sf_bessel_I1(x);
        double i2 = gsl_sf_bessel_In(2,x);

        return (exp(-gamma_relax*t)/(20.*M_PI))*(2.*gamma_relax*gamma_relax/c_diff)*(i1/(c_diff*u) + 4.*t*i2/u/u);
    }else{
        return 0.0;
    }
    
}

double CausalLiquefier::causal_diffusion_delta(double t, double r) const {

    double r_w = width_delta;
    if( c_diff*t <= width_delta ){
        r_w = c_diff*t;
    }
    
    if( r >= c_diff*t - r_w &&
        r < c_diff*t ){
        
        return (exp(-gamma_relax*t)/(20.*M_PI))*(8. - 3.*exp(-gamma_relax*t) + 2.*gamma_relax*t +4.*gamma_relax*gamma_relax*t*t )/r/r/r_w;
        
    }else{
        return 0.0;
    }
    
}

    
double CausalLiquefier::get_t(double tau, double eta)const{
    return tau*cosh(eta);
}
double CausalLiquefier::get_z(double tau, double eta)const{
        return tau*sinh(eta);
}
    
double CausalLiquefier::get_ptau(double pt, double pz, double eta)const{
        return pt*cosh(eta) - pz*cosh(eta);
}

double CausalLiquefier::get_peta(double pt, double pz, double eta)const{
        return pz*cosh(eta) - pt*cosh(eta);
}

    
    
};
