// Copyright @ Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <cmath>
#include <iostream>
#include "../src/fluid_dynamics.h"
#include "./brick_jetscape.h"

using namespace std;


Brick::Brick() {
    // initialize the parameter reader
    T_brick = 0.3;  // GeV

    hydro_status = NOT_START;
}


Brick::~Brick() {}


void Brick::initialize_hydro(Parameter parameter_list) {
    hydro_status = INITIALIZED;
}


void Brick::evolve_hydro() {
    hydro_status = FINISHED;
}


void Brick::get_hydro_info(real t, real x, real y, real z,
                           FluidCellInfo* fluid_cell_info_ptr) {
    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities
    fluid_cell_info_ptr->energy_density = 0.0;
    fluid_cell_info_ptr->entropy_density = 0.0;
    fluid_cell_info_ptr->temperature = T_brick;
    fluid_cell_info_ptr->pressure = 0.0;
    // QGP fraction
    fluid_cell_info_ptr->qgp_fraction = 1.0;
    // chemical potentials
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    // dynamical quantites
    fluid_cell_info_ptr->vx = 0.0;
    fluid_cell_info_ptr->vy = 0.0;
    fluid_cell_info_ptr->vz = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            fluid_cell_info_ptr->pi[i][j] = 0.0;
        }
    }
    fluid_cell_info_ptr->bulk_Pi = 0.0;
}
