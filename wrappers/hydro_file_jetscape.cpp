// Copyright @ Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <cmath>
#include <iostream>
#include "fluid_dynamics.h"
#include "hydro_file_jetscape.h"
#include "Hydroinfo_h5.h"
#include "Hydroinfo_MUSIC.h"
#include "ParameterReader.h"

using namespace std;


HydroFile::HydroFile() {
    // initialize the parameter reader
    paraRdr = new ParameterReader();
    paraRdr->readFromFile("parameters.dat");
    paraRdr->echo();

    T_c = 0.15;

    hydro_status = NOT_START;
}


HydroFile::~HydroFile() {
    if (hydro_status != NOT_START) {
        delete paraRdr;
    }
}


void HydroFile::initialize_hydro(Parameter parameter_list) {
    // this function loads the hydro files
    hydro_type = paraRdr->getVal("hydro_type");
    load_viscous = paraRdr->getVal("load_viscous_info");
    if (hydro_type == 1) {
        hydroinfo_h5_ptr = new HydroinfoH5("JetData.h5", 500, load_viscous);
    } else if (hydro_type == 2) {
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 8;
        int nskip_tau = paraRdr->getVal("hydro_nskip_tau");
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau);
    } else if (hydro_type == 3) {
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 9;
        int nskip_tau = paraRdr->getVal("hydro_nskip_tau");
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau);
    } else if (hydro_type == 4) {
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 10;
        int nskip_tau = 1;
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau);
    } else {
        cout << "main: unrecognized hydro_type = " << hydro_type << endl;
        exit(1);
    }
    hydro_status = INITIALIZED;
}


void HydroFile::evolve_hydro() {
    hydro_status = FINISHED;
}


void HydroFile::get_hydro_info(real t, real x, real y, real z,
                               FluidCellInfo* fluid_cell_info_ptr) {
    // this function returns the thermodynamic and dynamical information at
    // the given space-time point

    double t_local = static_cast<double>(t);
    double x_local = static_cast<double>(x);
    double y_local = static_cast<double>(y);
    double z_local = static_cast<double>(z);

    // initialize the 
    fluidCell *temp_fluid_cell_ptr = new fluidCell;
    if (hydro_type == 1) {  // for OSU 2+1d hydro
        double tau_local = sqrt(t*t - z*z);
        if (isnan(tau_local)) {  // check
            cout << "[Error]: HydroFile::get_hydro_info(): "
                 << "tau is nan!" << endl;
            cout << "please check: t = " << t << ", z = " << z << endl;
            exit(1);
        }
        hydroinfo_h5_ptr->getHydroinfo(tau_local, x_local, y_local,
                                       temp_fluid_cell_ptr);
    } else if (hydro_type == 2 || hydro_type == 3 || hydro_type == 4) {
        hydroinfo_MUSIC_ptr->getHydroValues(x_local, y_local, z_local, t_local,
                                            temp_fluid_cell_ptr);
    }

    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities
    fluid_cell_info_ptr->energy_density = (
                                static_cast<real>(temp_fluid_cell_ptr->ed));
    fluid_cell_info_ptr->entropy_density = (
                                static_cast<real>(temp_fluid_cell_ptr->sd));
    fluid_cell_info_ptr->temperature = (
                        static_cast<real>(temp_fluid_cell_ptr->temperature));
    fluid_cell_info_ptr->pressure = (
                        static_cast<real>(temp_fluid_cell_ptr->pressure));
    // QGP fraction
    double qgp_fraction_local = 1.0;
    if (temp_fluid_cell_ptr->temperature < T_c)
        qgp_fraction_local = 0.0;
    fluid_cell_info_ptr->qgp_fraction = static_cast<real>(qgp_fraction_local);
    // chemical potentials
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    // dynamical quantites
    fluid_cell_info_ptr->vx = static_cast<real>(temp_fluid_cell_ptr->vx);
    fluid_cell_info_ptr->vy = static_cast<real>(temp_fluid_cell_ptr->vy);
    fluid_cell_info_ptr->vz = static_cast<real>(temp_fluid_cell_ptr->vz);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            fluid_cell_info_ptr->pi[i][j] = (
                            static_cast<real>(temp_fluid_cell_ptr->pi[i][j]));
        }
    }
    fluid_cell_info_ptr->bulk_Pi = (
                            static_cast<real>(temp_fluid_cell_ptr->bulkPi));
    delete temp_fluid_cell_ptr;
}
