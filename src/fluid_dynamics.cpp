// Copyright JETSCAPE Collaboration @ 2016
// This is a general basic class for hydrodynamics
// This is written by Longgang Pang and Chun Shen

#include "./fluid_dynamics.h"


real FluidDynamics::get_energy_density(real time, real x, real y, real z) {
    // this function returns the energy density [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real energy_density = fluid_cell_ptr->energy_density;
    delete fluid_cell_ptr;
    return(energy_density);
}

real FluidDynamics::get_entropy_density(real time, real x, real y, real z) {
    // this function returns the entropy density [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real entropy_density = fluid_cell_ptr->entropy_density;
    delete fluid_cell_ptr;
    return(entropy_density);
}

real FluidDynamics::get_temperature(real time, real x, real y, real z) {
    // this function returns the temperature [GeV] at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real temperature = fluid_cell_ptr->temperature;
    delete fluid_cell_ptr;
    return(temperature);
}

real FluidDynamics::get_qgp_fraction(real time, real x, real y, real z) {
    // this function returns the QGP fraction at a space time point
    // (time, x, y, z)
    FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    get_hydro_info(time, x, y, z, fluid_cell_ptr);
    real qgp_fraction = fluid_cell_ptr->qgp_fraction;
    delete fluid_cell_ptr;
    return(qgp_fraction);
}
