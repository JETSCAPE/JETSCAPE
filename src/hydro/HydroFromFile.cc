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

#include <stdio.h>
#include <sys/stat.h>
#include <MakeUniqueHelper.h>

#include <cstring>
#include <sstream>
#include <cmath>
#include <iostream>

#include "JetScapeLogger.h"

#include "HydroFromFile.h"

using namespace Jetscape;

HydroFromFile::HydroFromFile() {
    hydro_status = NOT_START;
    SetId("hydroFromFile");
}


HydroFromFile::~HydroFromFile() {
    clean_hydro_event();
}


//! this function loads the hydro files
void HydroFromFile::initialize_hydro(Parameter parameter_list) {
    JSDEBUG << "Initialize hydro from file (Test) ...";
    VERBOSE(8);
    para_ = GetHydroXML()->FirstChildElement("hydro_from_file");
    if (!para_) {
        WARN << " : hydro_from_file not properly initialized in XML file ...";
        exit(-1);
    }

    string s = para_->FirstChildElement("name")->GetText();
    JSDEBUG << s << " to be initilizied ...";

    para_->FirstChildElement("hydro_type")->QueryIntText(&hydro_type_);
    para_->FirstChildElement("load_viscous_info")->QueryIntText(&load_viscous_);
    para_->FirstChildElement("read_hydro_every_ntau")->QueryIntText(
                                                                &nskip_tau_);
    para_->FirstChildElement("T_c")->QueryDoubleText(&T_c_);
    para_->FirstChildElement("read_in_multiple_hydro")->QueryIntText(
                                            &flag_read_in_multiple_hydro_);
    hydro_event_idx_ = 0;

    if (hydro_type_ == 1) {
#ifdef USE_HDF5
        hydroinfo_h5_ptr = new HydroinfoH5();
#else
        WARN << " : hydro_type == 1 requires the hdf5 library~";
        WARN << " : please check your inputs~";
        exit(-1);
#endif
    } else if (hydro_type_ == 2 || hydro_type_ == 3 || hydro_type_ == 4) {
        hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
    }

    hydro_status = INITIALIZED;
}

//! This function load a VISHNew hydro event
void HydroFromFile::read_in_hydro_event(string VISH_filename, int buffer_size,
                                    int load_viscous) {
    INFO << "read in a VISHNew hydro event from file " << VISH_filename;
    if (hydro_type_ == 1) {
#ifdef USE_HDF5
        hydroinfo_h5_ptr->readHydroinfoH5(VISH_filename, buffer_size,
                                          load_viscous);
#endif
    }
    hydro_status = FINISHED;
}


//! This function load a MUSIC hydro event
void HydroFromFile::read_in_hydro_event(string MUSIC_input_file,
                                    string MUSIC_hydro_ideal_file,
                                    int nskip_tau) {
    INFO << "read in a MUSIC hydro event from file " << MUSIC_hydro_ideal_file;
    if (hydro_type_ == 2) {
        int hydro_mode = 8;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    } else if (hydro_type_ == 3) {
        int hydro_mode = 9;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    } else if (hydro_type_ == 4) {
        int hydro_mode = 10;
        string hydro_shear_file = "";
        string hydro_bulk_file = "";
        hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau,
            MUSIC_input_file, MUSIC_hydro_ideal_file,
            hydro_shear_file, hydro_bulk_file);
    }
    hydro_status = FINISHED;
}


void HydroFromFile::evolve_hydro() {
    if (hydro_status == FINISHED) {
        clean_hydro_event();
        hydro_event_idx_ = ini->get_event_id();
    }

    if (hydro_type_ == 1) {
        string filename;
        if (flag_read_in_multiple_hydro_ == 0) {
            filename = para_->FirstChildElement("VISH_file")->GetText();
        } else {
            string folder = (
                    para_->FirstChildElement("hydro_files_folder")->GetText());
            std::ostringstream hydro_filename;
            hydro_filename << folder << "/event-" << hydro_event_idx_
                           << "/JetData.h5";
            filename = hydro_filename.str();
        }
#ifdef USE_HDF5
        read_in_hydro_event(filename, 500, load_viscous_);
#endif
        hydro_status = FINISHED;
    } else if (hydro_type_ == 2) {
        string input_file;
        string hydro_ideal_file;
        if (flag_read_in_multiple_hydro_ == 0) {
            input_file = (
                    para_->FirstChildElement("MUSIC_input_file")->GetText());
            hydro_ideal_file = (
                    para_->FirstChildElement("MUSIC_file")->GetText());
        } else {
            string folder = (
                    para_->FirstChildElement("hydro_files_folder")->GetText());
            std::ostringstream input_filename;
            std::ostringstream hydro_filename;
            input_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_input";
            hydro_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_evo.dat";
            input_file = input_filename.str();
            hydro_ideal_file = hydro_filename.str();
        }
        read_in_hydro_event(input_file, hydro_ideal_file, nskip_tau_);
    } else if (hydro_type_ == 3) {
        string input_file;
        string hydro_ideal_file;
        if (flag_read_in_multiple_hydro_ == 0) {
            input_file = (
                    para_->FirstChildElement("MUSIC_input_file")->GetText());
            hydro_ideal_file = (
                    para_->FirstChildElement("MUSIC_file")->GetText());
        } else {
            string folder = (
                    para_->FirstChildElement("hydro_files_folder")->GetText());
            std::ostringstream input_filename;
            std::ostringstream hydro_filename;
            input_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_input";
            hydro_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_evo.dat";
            input_file = input_filename.str();
            hydro_ideal_file = hydro_filename.str();
        }
        read_in_hydro_event(input_file, hydro_ideal_file, nskip_tau_);
    } else if (hydro_type_ == 4) {
        string input_file;
        string hydro_ideal_file;
        if (flag_read_in_multiple_hydro_ == 0) {
            input_file = (
                    para_->FirstChildElement("MUSIC_input_file")->GetText());
            hydro_ideal_file = (
                    para_->FirstChildElement("MUSIC_file")->GetText());
        } else {
            string folder = (
                    para_->FirstChildElement("hydro_files_folder")->GetText());
            std::ostringstream input_filename;
            std::ostringstream hydro_filename;
            input_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_input";
            hydro_filename << folder << "/event-" << hydro_event_idx_
                           << "/MUSIC_evo.dat";
            input_file = input_filename.str();
            hydro_ideal_file = hydro_filename.str();
        }
        read_in_hydro_event(input_file, hydro_ideal_file, 1);
    } else {
        WARN << "main: unrecognized hydro_type = " << hydro_type_;
        exit(1);
    }
}


//! clean up hydro event
void HydroFromFile::clean_hydro_event() {
    INFO << " clean up the loaded hydro event ...";
    if (hydro_type_ == 1) {
#ifdef USE_HDF5
        hydroinfo_h5_ptr->clean_hydro_event();
#endif
    } else {
        hydroinfo_MUSIC_ptr->clean_hydro_event();
    }
    hydro_status = NOT_START;
}


//! this function returns the thermodynamic and dynamical information at
//! the given space-time point
void HydroFromFile::get_hydro_info(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
        std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr) {
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
    if (hydro_type_ == 1) {  // for OSU 2+1d hydro
        double tau_local = sqrt(t*t - z*z);
        double eta_local = 0.5*log((t + z)/(t - z + 1e-15));
        if (std::isnan(tau_local)) {  // check
            WARN << "[Error]: HydroFromFile::get_hydro_info(): "
                 << "tau is nan!";
            WARN << "please check: t = " << t << ", z = " << z;
            exit(1);
        }
#ifdef USE_HDF5
        hydroinfo_h5_ptr->getHydroinfo(tau_local, x_local, y_local,
                                       temp_fluid_cell_ptr);
        // compute the flow velocity in the lab frame
        double u0_perp = (
            1./sqrt(1. - temp_fluid_cell_ptr->vx*temp_fluid_cell_ptr->vx
                       - temp_fluid_cell_ptr->vy*temp_fluid_cell_ptr->vy));
        double u0 = u0_perp*cosh(eta_local);
        temp_fluid_cell_ptr->vx *= u0_perp/u0;
        temp_fluid_cell_ptr->vy *= u0_perp/u0;
        temp_fluid_cell_ptr->vz = z/(t + 1e-15);
#endif
    } else if (hydro_type_ == 2 || hydro_type_ == 3 || hydro_type_ == 4) {
        hydroinfo_MUSIC_ptr->getHydroValues(x_local, y_local, z_local, t_local,
                                            temp_fluid_cell_ptr);
    }

    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities
    fluid_cell_info_ptr = std::make_unique<FluidCellInfo>();
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
    if (temp_fluid_cell_ptr->temperature < T_c_) {
        qgp_fraction_local = 0.0;
    }
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
