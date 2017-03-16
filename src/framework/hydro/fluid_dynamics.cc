// Copyright JETSCAPE Collaboration @ 2016
// This is a general basic class for hydrodynamics
// This is written by Longgang Pang and Chun Shen

#include <iostream>
#include "./fluid_dynamics.h"
#include "./linear_interpolation.h"

using namespace std;

/** For one given time step id_tau,
 * get FluidCellInfo at spatial point (x, y, eta)*/
FluidCellInfo EvolutionHistory::get_at_time_step(int id_tau,
                                     real x, real y, real eta) {
    int id_x = get_id_x(x);
    int id_y = get_id_y(y);
    int id_eta = get_id_eta(eta);
    // cijk for idx=i, idy=j and id_eta=k
    int c000 = cell_index(id_tau, id_x, id_y, id_eta);
    int c001 = cell_index(id_tau, id_x, id_y, id_eta+1);
    int c010 = cell_index(id_tau, id_x, id_y+1, id_eta);
    int c011 = cell_index(id_tau, id_x, id_y+1, id_eta+1);
    int c100 = cell_index(id_tau, id_x+1, id_y, id_eta);
    int c101 = cell_index(id_tau, id_x+1, id_y, id_eta+1);
    int c110 = cell_index(id_tau, id_x+1, id_y+1, id_eta);
    int c111 = cell_index(id_tau, id_x+1, id_y+1, id_eta+1);
    real x0 = x_coord(id_x);
    real x1 = x_coord(id_x + 1);
    real y0 = y_coord(id_y);
    real y1 = y_coord(id_y + 1);
    real eta0 = eta_coord(id_eta);
    real eta1 = eta_coord(id_eta + 1);

    return trilinear_int(x0, x1, y0, y1, eta0, eta1,
            data.at(c000), data.at(c001), data.at(c010), data.at(c011),
            data.at(c100), data.at(c101), data.at(c110), data.at(c111),
            x, y, eta);
}

// do interpolation along time direction; we may also need high order
// interpolation functions 
FluidCellInfo EvolutionHistory::get(real tau, real x, real y, real eta){
    check_in_range(tau, x, y, eta);
    int id_tau = get_id_tau(tau);
    real tau0 = tau_coord(id_tau);
    real tau1 = tau_coord(id_tau + 1);
    FluidCellInfo bulk0 = get_at_time_step(id_tau, x, y, eta);
    FluidCellInfo bulk1 = get_at_time_step(id_tau+1, x, y, eta);
    return linear_int(tau0, tau1, bulk0, bulk1, tau);
}


real FluidDynamicsBase::get_energy_density(real time, real x, real y, real z) {
    // this function returns the energy density [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real energy_density = fluid_cell_ptr->energy_density;
    delete fluid_cell_ptr;
    return(energy_density);
}

real FluidDynamicsBase::get_entropy_density(real time, real x, real y, real z) {
    // this function returns the entropy density [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real entropy_density = fluid_cell_ptr->entropy_density;
    delete fluid_cell_ptr;
    return(entropy_density);
}

real FluidDynamicsBase::get_temperature(real time, real x, real y, real z) {
    // this function returns the temperature [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real temperature = fluid_cell_ptr->temperature;
    delete fluid_cell_ptr;
    return(temperature);
}

real FluidDynamicsBase::get_qgp_fraction(real time, real x, real y, real z) {
    // this function returns the QGP fraction at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real qgp_fraction = fluid_cell_ptr->qgp_fraction;
    delete fluid_cell_ptr;
    return(qgp_fraction);
}

void FluidDynamicsBase::print_fluid_cell_information(
                                        FluidCellInfo* fluid_cell_info_ptr) {
    // this function print out the information of the fluid cell to the screen
    cout << "=======================================================" << endl;
    cout << "print out cell information:" << endl;
    cout << "=======================================================" << endl;
    cout << "energy density = " << fluid_cell_info_ptr->energy_density
         << " GeV/fm^3." << endl;
    cout << "entropy density = " << fluid_cell_info_ptr->entropy_density
         << " 1/fm^3." << endl;
    cout << "temperature = " << fluid_cell_info_ptr->temperature
         << " GeV." << endl;
    cout << "pressure = " << fluid_cell_info_ptr->pressure
         << " GeV/fm^3." << endl;
    cout << "QGP_fraction = " << fluid_cell_info_ptr->qgp_fraction
         << endl;
    cout << "mu_B = " << fluid_cell_info_ptr->mu_B << " GeV." << endl;
    cout << "mu_S = " << fluid_cell_info_ptr->mu_S << " GeV." << endl;
    cout << "mu_C = " << fluid_cell_info_ptr->mu_C << " GeV." << endl;
    cout << "vx = " << fluid_cell_info_ptr->vx << endl;
    cout << "vy = " << fluid_cell_info_ptr->vy << endl;
    cout << "vz = " << fluid_cell_info_ptr->vz << endl;
    cout << "shear viscous pi^{munu} (GeV/fm^3): " << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << fluid_cell_info_ptr->pi[i][j] << "   ";
        }
        cout << endl;
    }
    cout << "bulk_Pi = " << fluid_cell_info_ptr->bulk_Pi << " GeV/fm^3"
         << endl;
    cout << "=======================================================" << endl;
}

void FluidCellInfo::Print()
{
    // this function print out the information of the fluid cell to the screen
    cout << "=======================================================" << endl;
    cout << "print out cell information:" << endl;
    cout << "=======================================================" << endl;
    cout << "energy density = " << energy_density
         << " GeV/fm^3." << endl;
    cout << "entropy density = " << entropy_density
         << " 1/fm^3." << endl;
    cout << "temperature = " << temperature
         << " GeV." << endl;
    cout << "pressure = " << pressure
         << " GeV/fm^3." << endl;
    cout << "QGP_fraction = " << qgp_fraction
         << endl;
    cout << "mu_B = " << mu_B << " GeV." << endl;
    cout << "mu_S = " << mu_S << " GeV." << endl;
    cout << "mu_C = " << mu_C << " GeV." << endl;
    cout << "vx = " << vx << endl;
    cout << "vy = " << vy << endl;
    cout << "vz = " << vz << endl;
    cout << "shear viscous pi^{munu} (GeV/fm^3): " << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << pi[i][j] << "   ";
        }
        cout << endl;
    }
    cout << "bulk_Pi = " << bulk_Pi << " GeV/fm^3"
         << endl;
    cout << "=======================================================" << endl;
}
