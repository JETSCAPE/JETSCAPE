/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// Copyright @ Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <cmath>
#include <iostream>
#include <MakeUniqueHelper.h>
#include "GubserHydro.h"

#define hbarc 0.19733

// using namespace std;
using namespace Jetscape;

GubserHydro::GubserHydro() : FluidDynamics(){
    // initialize the parameter reader
    q = 1.0;
    e_0 = 1.0;

    hydro_status = NOT_START;
    SetId("Gubser");
    VERBOSE(8);
}


GubserHydro::~GubserHydro() { VERBOSE(8);}


void GubserHydro::initialize_hydro(Parameter parameter_list) {
   VERBOSE(8);
    hydro_status = INITIALIZED;
}


void GubserHydro::evolve_hydro() {
   VERBOSE(8);
    hydro_status = FINISHED;
}


double GubserHydro::temperature(double e_local) {
    double N_c = 3.;
    double N_f = 2.5;
    double T_local = pow(
        90.0/M_PI/M_PI*e_local/3./(2.*(N_c*N_c - 1.) + 7./2.*N_c*N_f), 0.25);
    T_local *= hbarc;
    return(T_local);
}


void GubserHydro::get_hydro_info(real t, real x, real y, real z,
//                                  FluidCellInfo* fluid_cell_info_ptr) {
				 std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr){
    // create the unique FluidCellInfo here
    fluid_cell_info_ptr=std::make_unique<FluidCellInfo>();
    
    double t_local = static_cast<double>(t);
    double x_local = static_cast<double>(x);
    double y_local = static_cast<double>(y);
    double z_local = static_cast<double>(z);
        
    double tau_local = sqrt(t*t - z*z);
    double r_local = sqrt(x_local*x_local + y_local*y_local);

    double temp = (1. + 2.*q*q*(tau_local*tau_local + r_local*r_local)
                   + q*q*q*q*pow(tau_local*tau_local - r_local*r_local, 2));

    double e_local = (
        (e_0/pow(tau_local, 4./3.))*(pow(2.*q, 8./3.))/(pow(temp, 4./3.)));
    double T_local = temperature(e_local);           // GeV
    e_local *= hbarc;                                // GeV/fm^3
    double p_local = e_local/3.;                     // GeV/fm^3
    double s_local = (e_local + p_local)/T_local;    // 1/fm^3

    double kappa = atanh(
        (2.*q*q*tau_local*r_local)/(1. + q*q*tau_local*tau_local
                                    + q*q*r_local*r_local));
    double ux_local = sinh(kappa)*x_local/(r_local + 1e-15);
    double uy_local = sinh(kappa)*y_local/(r_local + 1e-15);
    double gamma = sqrt(1. + ux_local*ux_local + uy_local*uy_local);
    double vx_local = ux_local/gamma;
    double vy_local = uy_local/gamma;
    double vz_local = z/t;

    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities
    fluid_cell_info_ptr->energy_density = e_local;
    fluid_cell_info_ptr->entropy_density = s_local;
    fluid_cell_info_ptr->temperature = T_local;
    fluid_cell_info_ptr->pressure = p_local;
    // QGP fraction
    fluid_cell_info_ptr->qgp_fraction = 1.0;
    // chemical potentials
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    // dynamical quantites
    fluid_cell_info_ptr->vx = vx_local;
    fluid_cell_info_ptr->vy = vy_local;
    fluid_cell_info_ptr->vz = vz_local;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            fluid_cell_info_ptr->pi[i][j] = 0.0;
        }
    }
    fluid_cell_info_ptr->bulk_Pi = 0.0;
}
