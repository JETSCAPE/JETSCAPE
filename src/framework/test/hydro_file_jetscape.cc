// Copyright @ Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <cmath>
#include <iostream>

#include "JetScapeLogger.h"

#include "hydro_file_jetscape.h"

using namespace Jetscape;

HydroFile::HydroFile() {
    hydro_status = NOT_START;
    SetId("hydroFromFile");
}


HydroFile::~HydroFile() {
    clean_hydro_event();
}


//! this function loads the hydro files
void HydroFile::initialize_hydro(Parameter parameter_list) {
    JSDEBUG << "Initialize hydro from file (Test) ...";
    VERBOSE(8);
    tinyxml2::XMLElement *para =
                    GetHydroXML()->FirstChildElement("hydro_from_file");
    if (!para) {
        WARN << " : hydro_from_file not properly initialized in XML file ...";
        exit(-1);
    }

    string s = para->FirstChildElement("name")->GetText();
    JSDEBUG << s << " to be initilizied ...";

    para->FirstChildElement("hydro_type")->QueryIntText(&hydro_type);
    para->FirstChildElement("load_viscous_info")->QueryBoolText(&load_viscous);
    para->FirstChildElement("T_c")->QueryDoubleText(&T_c);
    hydro_status = INITIALIZED;

    int nskip_tau = 1;
    if (hydro_type == 1) {
        string filename = para->FirstChildElement("VISH_file")->GetText();
#ifdef USE_HDF5
        hydroinfo_h5_ptr = new HydroinfoH5(filename, 500, load_viscous);
#else
        WARN << " : hydro_type == 1 requires the hdf5 library~";
        WARN << " : please check your inputs~";
        exit(-1);
#endif
    } else if (hydro_type == 2) {
        string input_file =
                para->FirstChildElement("MUSIC_input_file")->GetText();
        string hydro_ideal_file =
                para->FirstChildElement("MUSIC_file")->GetText();
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        para->FirstChildElement("read_hydro_every_ntau")->QueryIntText(
                                                                &nskip_tau);
        read_in_hydro_event(input_file, hydro_ideal_file, nskip_tau);
    } else if (hydro_type == 3) {
        string input_file =
                para->FirstChildElement("MUSIC_input_file")->GetText();
        string hydro_ideal_file =
                para->FirstChildElement("MUSIC_file")->GetText();
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        para->FirstChildElement("read_hydro_every_ntau")->QueryIntText(
                                                                &nskip_tau);
        read_in_hydro_event(input_file, hydro_ideal_file, nskip_tau);
    } else if (hydro_type == 4) {
        string input_file =
                para->FirstChildElement("MUSIC_input_file")->GetText();
        string hydro_ideal_file =
                para->FirstChildElement("MUSIC_file")->GetText();
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
        nskip_tau = 1;
        read_in_hydro_event(input_file, hydro_ideal_file, nskip_tau);
    } else {
        WARN << "main: unrecognized hydro_type = " << hydro_type;
        exit(1);
    }
}

     
//! This function load a VISHNew hydro event
void HydroFile::read_in_hydro_event(string VISH_filename, int buffer_size,
                                    int load_viscous) {
    INFO << "read in a VISHNew hydro event from file " << VISH_filename;
    if (hydro_type == 1) {
#ifdef USE_HDF5
        hydroinfo_h5_ptr->readHydroinfoH5(VISH_filename, buffer_size,
                                          load_viscous);
#else
        WARN << " : hydro_type == 1 requires the hdf5 library~";
        WARN << " : please check your inputs~";
        exit(-1);
#endif
    }
    hydro_status = FINISHED;
}


//! This function load a MUSIC hydro event
void HydroFile::read_in_hydro_event(string MUSIC_input_file,
                                    string MUSIC_hydro_ideal_file,
                                    int nskip_tau) {
    INFO << "read in a MUSIC hydro event from file " << MUSIC_hydro_ideal_file;
    if (hydro_type == 2) {
        int hydro_mode = 8;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    } else if (hydro_type == 3) {
        int hydro_mode = 9;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    } else if (hydro_type == 4) {
        int hydro_mode = 10;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    }
    hydro_status = FINISHED;
}


void HydroFile::evolve_hydro() {
    VERBOSE(8);
    hydro_status = FINISHED;
}


//! clean up hydro event
void HydroFile::clean_hydro_event() {
    INFO << " clean up the loaded hydro event ...";
    if (hydro_type == 1) {
#ifdef USE_HDF5
        hydroinfo_h5_ptr->clean_hydro_event();
#else
        WARN << " : hydro_type == 1 requires the hdf5 library~";
        WARN << " : please check your inputs~";
        exit(-1);
#endif
    } else {
        hydroinfo_MUSIC_ptr->clean_hydro_event();
    }
    hydro_status = NOT_START;
}


//! this function returns the thermodynamic and dynamical information at
//! the given space-time point
void HydroFile::get_hydro_info(real t, real x, real y, real z,
                               FluidCellInfo* fluid_cell_info_ptr) {
    if (hydro_status != FINISHED) {
        WARN << "Hydro not run yet ...";
        exit(-1);
    }

    double t_local = static_cast<double>(t);
    double x_local = static_cast<double>(x);
    double y_local = static_cast<double>(y);
    double z_local = static_cast<double>(z);

    // initialize the fluid cell pointer
    fluidCell *temp_fluid_cell_ptr = new fluidCell;
    if (hydro_type == 1) {  // for OSU 2+1d hydro
        double tau_local = sqrt(t*t - z*z);
        if (std::isnan(tau_local)) {  // check
            WARN << "[Error]: HydroFile::get_hydro_info(): "
                 << "tau is nan!";
            WARN << "please check: t = " << t << ", z = " << z;
            exit(1);
        }
#ifdef USE_HDF5
        hydroinfo_h5_ptr->getHydroinfo(tau_local, x_local, y_local,
                                       temp_fluid_cell_ptr);
#endif
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
