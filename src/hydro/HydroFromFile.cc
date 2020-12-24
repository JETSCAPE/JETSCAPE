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

// Register the module with the base class
RegisterJetScapeModule<HydroFromFile> HydroFromFile::reg("HydroFromFile");

HydroFromFile::HydroFromFile() {
  hydro_status = NOT_START;
  SetId("hydroFromFile");
}

HydroFromFile::~HydroFromFile() { clean_hydro_event(); }

//! this function loads the hydro files
void HydroFromFile::InitializeHydro(Parameter parameter_list) {
  JSDEBUG << "Initialize hydro from file (Test) ...";
  VERBOSE(8);

  string s = GetXMLElementText({"Hydro", "hydro_from_file", "name"});
  JSDEBUG << s << " to be initilizied ...";

  hydro_type_ = GetXMLElementInt({"Hydro", "hydro_from_file", "hydro_type"});
  load_viscous_ =
      GetXMLElementInt({"Hydro", "hydro_from_file", "load_viscous_info"});
  nskip_tau_ =
      GetXMLElementInt({"Hydro", "hydro_from_file", "read_hydro_every_ntau"});
  T_c_ = GetXMLElementDouble({"Hydro", "hydro_from_file", "T_c"});
  flag_read_in_multiple_hydro_ =
      GetXMLElementInt({"Hydro", "hydro_from_file", "read_in_multiple_hydro"});

  hydro_event_idx_ = 0;

  if (hydro_type_ == 1) {
#ifdef USE_HDF5
    hydroinfo_h5_ptr = new HydroinfoH5();
#else
    JSWARN << " : hydro_type == 1 requires the hdf5 library~";
    JSWARN << " : please check your inputs~";
    exit(-1);
#endif
  } else if (hydro_type_ == 2 || hydro_type_ == 3 || hydro_type_ == 4) {
    hydroinfo_MUSIC_ptr = new Hydroinfo_MUSIC();
    int verbose = GetXMLElementInt({"vlevel"});
    hydroinfo_MUSIC_ptr->set_verbose(verbose);
  }

  hydro_status = INITIALIZED;
}

//! This function load a VISHNew hydro event
void HydroFromFile::read_in_hydro_event(string VISH_filename, int buffer_size,
                                        int load_viscous) {
  JSINFO << "read in a VISHNew hydro event from file " << VISH_filename;
  if (hydro_type_ == 1) {
#ifdef USE_HDF5
    hydroinfo_h5_ptr->readHydroinfoH5(VISH_filename, buffer_size, load_viscous);
#endif
  }
  hydro_status = FINISHED;
}

//! This function load a MUSIC hydro event
void HydroFromFile::read_in_hydro_event(string MUSIC_input_file,
                                        string MUSIC_hydro_ideal_file,
                                        int nskip_tau) {
  JSINFO << "read in a MUSIC hydro event from file " << MUSIC_hydro_ideal_file;
  if (hydro_type_ == 2) {
    int hydro_mode = 8;
    string hydro_shear_file = "";
    string hydro_bulk_file = "";
    hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau, MUSIC_input_file,
                                       MUSIC_hydro_ideal_file, hydro_shear_file,
                                       hydro_bulk_file);
  } else if (hydro_type_ == 3) {
    int hydro_mode = 9;
    string hydro_shear_file = "";
    string hydro_bulk_file = "";
    hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau, MUSIC_input_file,
                                       MUSIC_hydro_ideal_file, hydro_shear_file,
                                       hydro_bulk_file);
  } else if (hydro_type_ == 4) {
    int hydro_mode = 10;
    string hydro_shear_file = "";
    string hydro_bulk_file = "";
    hydroinfo_MUSIC_ptr->readHydroData(hydro_mode, nskip_tau, MUSIC_input_file,
                                       MUSIC_hydro_ideal_file, hydro_shear_file,
                                       hydro_bulk_file);
  }
  hydro_status = FINISHED;
}

void HydroFromFile::EvolveHydro() {
  if (hydro_status == FINISHED) {
    clean_hydro_event();
    hydro_event_idx_ = ini->GetEventId();
  }

  if (hydro_type_ == 1) {
    string filename;
    if (flag_read_in_multiple_hydro_ == 0) {
      filename = GetXMLElementText({"Hydro", "hydro_from_file", "VISH_file"});
    } else {
      string folder =
          GetXMLElementText({"Hydro", "hydro_from_file", "hydro_files_folder"});
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
      input_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_input_file"});
      hydro_ideal_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_file"});
    } else {
      string folder = GetXMLElementText(
              {"Hydro", "hydro_from_file", "hydro_files_folder"});
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
      input_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_input_file"});
      hydro_ideal_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_file"});
    } else {
      string folder = GetXMLElementText(
              {"Hydro", "hydro_from_file", "hydro_files_folder"});
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
      input_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_input_file"});
      hydro_ideal_file = GetXMLElementText(
              {"Hydro", "hydro_from_file", "MUSIC_file"});
    } else {
      string folder = GetXMLElementText(
              {"Hydro", "hydro_from_file", "hydro_files_folder"});
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
    JSWARN << "main: unrecognized hydro_type = " << hydro_type_;
    exit(1);
  }
}

//! clean up hydro event
void HydroFromFile::clean_hydro_event() {
  JSINFO << " clean up the loaded hydro event ...";
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
void HydroFromFile::GetHydroInfo(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  if (hydro_status != FINISHED) {
    JSWARN << "Hydro not run yet ...";
    exit(-1);
  }

  double t_local = static_cast<double>(t);
  double x_local = static_cast<double>(x);
  double y_local = static_cast<double>(y);
  double z_local = static_cast<double>(z);

  // initialize the fluid cell pointer
  hydrofluidCell *temp_fluid_cell_ptr = new hydrofluidCell;
  if (hydro_type_ == 1) { // for OSU 2+1d hydro
    double tau_local = sqrt(t * t - z * z);
    double eta_local = 0.5 * log((t + z) / (t - z + 1e-15));
    if (std::isnan(tau_local)) { // check
      JSWARN << "[Error]: HydroFromFile::GetHydroInfo(): "
             << "tau is nan!";
      JSWARN << "please check: t = " << t << ", z = " << z;
      exit(1);
    }
#ifdef USE_HDF5
    hydroinfo_h5_ptr->getHydroinfo(tau_local, x_local, y_local,
                                   temp_fluid_cell_ptr);
    // compute the flow velocity in the lab frame
    double u0_perp =
        (1. / sqrt(1. - temp_fluid_cell_ptr->vx * temp_fluid_cell_ptr->vx -
                   temp_fluid_cell_ptr->vy * temp_fluid_cell_ptr->vy));
    double u0 = u0_perp * cosh(eta_local);
    temp_fluid_cell_ptr->vx *= u0_perp / u0;
    temp_fluid_cell_ptr->vy *= u0_perp / u0;
    temp_fluid_cell_ptr->vz = z / (t + 1e-15);
#endif
  } else if (hydro_type_ == 2 || hydro_type_ == 3 || hydro_type_ == 4) {
    hydroinfo_MUSIC_ptr->getHydroValues(x_local, y_local, z_local, t_local,
                                        temp_fluid_cell_ptr);
  }

  // assign all the quantites to JETSCAPE output
  // thermodyanmic quantities
  fluid_cell_info_ptr = make_unique<FluidCellInfo>();
  fluid_cell_info_ptr->energy_density =
      (static_cast<Jetscape::real>(temp_fluid_cell_ptr->ed));
  fluid_cell_info_ptr->entropy_density =
      (static_cast<Jetscape::real>(temp_fluid_cell_ptr->sd));
  fluid_cell_info_ptr->temperature =
      (static_cast<Jetscape::real>(temp_fluid_cell_ptr->temperature));
  fluid_cell_info_ptr->pressure =
      (static_cast<Jetscape::real>(temp_fluid_cell_ptr->pressure));
  // QGP fraction
  double qgp_fraction_local = 1.0;
  if (temp_fluid_cell_ptr->temperature < T_c_) {
    qgp_fraction_local = 0.0;
  }
  fluid_cell_info_ptr->qgp_fraction =
      static_cast<Jetscape::real>(qgp_fraction_local);
  // chemical potentials
  fluid_cell_info_ptr->mu_B = 0.0;
  fluid_cell_info_ptr->mu_C = 0.0;
  fluid_cell_info_ptr->mu_S = 0.0;
  // dynamical quantites
  fluid_cell_info_ptr->vx =
      static_cast<Jetscape::real>(temp_fluid_cell_ptr->vx);
  fluid_cell_info_ptr->vy =
      static_cast<Jetscape::real>(temp_fluid_cell_ptr->vy);
  fluid_cell_info_ptr->vz =
      static_cast<Jetscape::real>(temp_fluid_cell_ptr->vz);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      fluid_cell_info_ptr->pi[i][j] =
          (static_cast<Jetscape::real>(temp_fluid_cell_ptr->pi[i][j]));
    }
  }
  fluid_cell_info_ptr->bulk_Pi =
      (static_cast<Jetscape::real>(temp_fluid_cell_ptr->bulkPi));
  delete temp_fluid_cell_ptr;
}

double HydroFromFile::GetEventPlaneAngle() {
  double v2 = 0.0;
  double psi_2 = 0.0;
  if (hydro_type_ == 1) {
    std::ostringstream angle_filename;
    string folder =
        GetXMLElementText({"Hydro", "hydro_from_file", "hydro_files_folder"});
    angle_filename << folder << "/event-" << hydro_event_idx_
                   << "/EventPlanesFrzout.dat";
    std::ifstream inputfile(angle_filename.str().c_str());
    string dummy;
    std::getline(inputfile, dummy);
    std::getline(inputfile, dummy);
    std::getline(inputfile, dummy);
    inputfile >> dummy >> v2 >> dummy >> psi_2;
    inputfile.close();
  }
  return (psi_2);
}
