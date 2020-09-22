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

#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "JetScapeLogger.h"
#include "MusicWrapper.h"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<MpiMusic> MpiMusic::reg("MUSIC");

MpiMusic::MpiMusic() {
  hydro_status = NOT_START;
  freezeout_temperature = 0.0;
  doCooperFrye = 0;
  flag_output_evo_to_file = 0;
  has_source_terms = false;
  SetId("MUSIC");
  hydro_source_terms_ptr =
      std::shared_ptr<HydroSourceJETSCAPE>(new HydroSourceJETSCAPE());
}

MpiMusic::~MpiMusic() {}

void MpiMusic::InitializeHydro(Parameter parameter_list) {
  JSINFO << "Initialize MUSIC ...";
  VERBOSE(8);

  string input_file = GetXMLElementText({"Hydro", "MUSIC", "MUSIC_input_file"});
  doCooperFrye =
      GetXMLElementInt({"Hydro", "MUSIC", "Perform_CooperFrye_Feezeout"});

  music_hydro_ptr = std::unique_ptr<MUSIC>(new MUSIC(input_file));

  // overwrite input options
  flag_output_evo_to_file = (
      GetXMLElementInt({"Hydro", "MUSIC", "output_evolution_to_file"}));
  music_hydro_ptr->set_parameter("output_movie_flag",
                                 static_cast<double>(flag_output_evo_to_file));
  double tau_hydro = (
          GetXMLElementDouble({"Hydro", "MUSIC", "Initial_time_tau_0"}));
  music_hydro_ptr->set_parameter("Initial_time_tau_0", tau_hydro);

  double eta_over_s =
      GetXMLElementDouble({"Hydro", "MUSIC", "shear_viscosity_eta_over_s"});
  if (eta_over_s > 1e-6) {
    music_hydro_ptr->set_parameter("Viscosity_Flag_Yes_1_No_0", 1);
    music_hydro_ptr->set_parameter("Include_Shear_Visc_Yes_1_No_0", 1);
    music_hydro_ptr->set_parameter("Shear_to_S_ratio", eta_over_s);
  } else if (eta_over_s >= 0.) {
    music_hydro_ptr->set_parameter("Viscosity_Flag_Yes_1_No_0", 0);
    music_hydro_ptr->set_parameter("Include_Shear_Visc_Yes_1_No_0", 0);
  } else {
    JSWARN << "The input shear viscosity is negative! eta/s = " << eta_over_s;
    exit(1);
  }

  int flag_shear_Tdep = (
      GetXMLElementInt({"Hydro", "MUSIC", "T_dependent_Shear_to_S_ratio"}));
  if (flag_shear_Tdep > 0) {
    music_hydro_ptr->set_parameter("Viscosity_Flag_Yes_1_No_0", 1);
    if (flag_shear_Tdep == 3) {
      double shear_kinkT = (
        GetXMLElementDouble({"Hydro", "MUSIC", "eta_over_s_T_kink_in_GeV"}));
      music_hydro_ptr->set_parameter("eta_over_s_T_kink_in_GeV", shear_kinkT);
      double shear_lowTslope = (
        GetXMLElementDouble({"Hydro", "MUSIC",
                             "eta_over_s_low_T_slope_in_GeV"}));
      music_hydro_ptr->set_parameter("eta_over_s_low_T_slope_in_GeV",
                                     shear_lowTslope);
      double shear_highTslope = (
        GetXMLElementDouble({"Hydro", "MUSIC",
                             "eta_over_s_high_T_slope_in_GeV"}));
      music_hydro_ptr->set_parameter("eta_over_s_high_T_slope_in_GeV",
                                     shear_highTslope);
      double shear_kink = (
        GetXMLElementDouble({"Hydro", "MUSIC", "eta_over_s_at_kink"}));
      music_hydro_ptr->set_parameter("eta_over_s_at_kink", shear_kink);
    }
  }

  int flag_bulkvis = GetXMLElementInt(
          {"Hydro", "MUSIC", "temperature_dependent_bulk_viscosity"});
  if (flag_bulkvis != 0) {
    music_hydro_ptr->set_parameter("Include_Bulk_Visc_Yes_1_No_0", 1);
    music_hydro_ptr->set_parameter("T_dependent_Bulk_to_S_ratio",
                                   flag_bulkvis);
    if (flag_bulkvis == 3) {
        double bulk_max = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_max"});
        music_hydro_ptr->set_parameter("zeta_over_s_max", bulk_max);
        double bulk_peakT = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_T_peak_in_GeV"});
        music_hydro_ptr->set_parameter("zeta_over_s_T_peak_in_GeV",
                                       bulk_peakT);
        double bulk_width = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_width_in_GeV"});
        music_hydro_ptr->set_parameter("zeta_over_s_width_in_GeV", bulk_width);
        double bulk_asy = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_lambda_asymm"});
        music_hydro_ptr->set_parameter("zeta_over_s_lambda_asymm", bulk_asy);
    }
  }

  int flag_secondorderTerms = GetXMLElementInt(
          {"Hydro", "MUSIC", "Include_second_order_terms"});
  if (flag_secondorderTerms == 1) {
    music_hydro_ptr->set_parameter("Include_second_order_terms", 1);
  }

  freezeout_temperature =
      GetXMLElementDouble({"Hydro", "MUSIC", "freezeout_temperature"});
  if (freezeout_temperature > 0.05) {
    music_hydro_ptr->set_parameter("T_freeze", freezeout_temperature);
  } else {
    JSWARN << "The input freeze-out temperature is too low! T_frez = "
           << freezeout_temperature << " GeV!";
    exit(1);
  }


  music_hydro_ptr->add_hydro_source_terms(hydro_source_terms_ptr);
}

void MpiMusic::EvolveHydro() {
  VERBOSE(8);
  JSINFO << "Initialize density profiles in MUSIC ...";
  std::vector<double> entropy_density = ini->GetEntropyDensityDistribution();
  double dx = ini->GetXStep();
  double dz = ini->GetZStep();
  double z_max = ini->GetZMax();
  int nz = ini->GetZSize();
  if (pre_eq_ptr == nullptr) {
    JSWARN << "Missing the pre-equilibrium module ...";
  } else {
    music_hydro_ptr->initialize_hydro_from_jetscape_preequilibrium_vectors(
        dx, dz, z_max, nz, pre_eq_ptr->e_, pre_eq_ptr->utau_, pre_eq_ptr->ux_,
        pre_eq_ptr->uy_, pre_eq_ptr->ueta_, pre_eq_ptr->pi00_,
        pre_eq_ptr->pi01_, pre_eq_ptr->pi02_, pre_eq_ptr->pi03_,
        pre_eq_ptr->pi11_, pre_eq_ptr->pi12_, pre_eq_ptr->pi13_,
        pre_eq_ptr->pi22_, pre_eq_ptr->pi23_, pre_eq_ptr->pi33_,
        pre_eq_ptr->bulk_Pi_);
  }

  JSINFO << "initial density profile dx = " << dx << " fm";
  hydro_status = INITIALIZED;
  JSINFO << "number of source terms: "
         << hydro_source_terms_ptr->get_number_of_sources()
         << ", total E = " << hydro_source_terms_ptr->get_total_E_of_sources()
         << " GeV.";

  has_source_terms = false;
  if (hydro_source_terms_ptr->get_number_of_sources() > 0) {
    has_source_terms = true;
  }

  if (hydro_status == INITIALIZED) {
    JSINFO << "running MUSIC ...";
    music_hydro_ptr->run_hydro();
    hydro_status = FINISHED;
  }

  if (flag_output_evo_to_file == 1) {
    if (!has_source_terms) {
      // only the first hydro without source term will be stored
      // in memory for jet energy loss calculations
      PassHydroEvolutionHistoryToFramework();
      JSINFO << "number of fluid cells received by the JETSCAPE: "
             << bulk_info.data.size();
    }
    music_hydro_ptr->clear_hydro_info_from_memory();

    // add hydro_id to the hydro evolution filename
    std::ostringstream system_command;
    system_command << "mv evolution_for_movie_xyeta.dat "
                   << "evolution_for_movie_xyeta_" << GetId() << ".dat";
    system(system_command.str().c_str());

    //std::vector<SurfaceCellInfo> surface_cells;
    //if (freezeout_temperature > 0.0) {
    //  FindAConstantTemperatureSurface(freezeout_temperature, surface_cells);
    //}
  }

  collect_freeze_out_surface();

  if (hydro_status == FINISHED && doCooperFrye == 1) {
    music_hydro_ptr->run_Cooper_Frye();
  }
}

void MpiMusic::collect_freeze_out_surface() {
  system("rm surface.dat 2> /dev/null");

  std::ostringstream surface_filename;
  surface_filename << "surface_" << GetId() << ".dat";

  std::ostringstream system_command;
  system_command << "rm " << surface_filename.str() << " 2> /dev/null";
  system(system_command.str().c_str());
  system_command.str("");
  system_command.clear();
  system_command << "cat surface_eps* >> " << surface_filename.str();
  system(system_command.str().c_str());
  system_command.str("");
  system_command.clear();

  system_command << "ln -s " << surface_filename.str() << " surface.dat";
  system(system_command.str().c_str());
  system_command.str("");
  system_command.clear();
  system("rm surface_eps* 2> /dev/null");
}

void MpiMusic::SetHydroGridInfo() {
  bulk_info.neta = music_hydro_ptr->get_neta();
  bulk_info.ntau = music_hydro_ptr->get_ntau();
  bulk_info.nx = music_hydro_ptr->get_nx();
  bulk_info.ny = music_hydro_ptr->get_nx();
  bulk_info.tau_min = music_hydro_ptr->get_hydro_tau0();
  bulk_info.dtau = music_hydro_ptr->get_hydro_dtau();
  bulk_info.x_min = -music_hydro_ptr->get_hydro_x_max();
  bulk_info.dx = music_hydro_ptr->get_hydro_dx();
  bulk_info.y_min = -music_hydro_ptr->get_hydro_x_max();
  bulk_info.dy = music_hydro_ptr->get_hydro_dx();
  bulk_info.eta_min = -music_hydro_ptr->get_hydro_eta_max();
  bulk_info.deta = music_hydro_ptr->get_hydro_deta();

  bulk_info.boost_invariant = music_hydro_ptr->is_boost_invariant();
}

void MpiMusic::PassHydroEvolutionHistoryToFramework() {
  clear_up_evolution_data();

  JSINFO << "Passing hydro evolution information to JETSCAPE ... ";
  auto number_of_cells = music_hydro_ptr->get_number_of_fluid_cells();
  JSINFO << "total number of fluid cells: " << number_of_cells;

  SetHydroGridInfo();

  fluidCell *fluidCell_ptr = new fluidCell;
  for (int i = 0; i < number_of_cells; i++) {
    std::unique_ptr<FluidCellInfo> fluid_cell_info_ptr(new FluidCellInfo);
    music_hydro_ptr->get_fluid_cell_with_index(i, fluidCell_ptr);

    fluid_cell_info_ptr->energy_density = fluidCell_ptr->ed;
    fluid_cell_info_ptr->entropy_density = fluidCell_ptr->sd;
    fluid_cell_info_ptr->temperature = fluidCell_ptr->temperature;
    fluid_cell_info_ptr->pressure = fluidCell_ptr->pressure;
    fluid_cell_info_ptr->vx = fluidCell_ptr->vx;
    fluid_cell_info_ptr->vy = fluidCell_ptr->vy;
    fluid_cell_info_ptr->vz = fluidCell_ptr->vz;
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    fluid_cell_info_ptr->qgp_fraction = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        fluid_cell_info_ptr->pi[i][j] = fluidCell_ptr->pi[i][j];
      }
    }
    fluid_cell_info_ptr->bulk_Pi = fluidCell_ptr->bulkPi;
    StoreHydroEvolutionHistory(fluid_cell_info_ptr);
  }
  delete fluidCell_ptr;
}

void MpiMusic::GetHydroInfo(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  GetHydroInfo_JETSCAPE(t, x, y, z, fluid_cell_info_ptr);
  //GetHydroInfo_MUSIC(t, x, y, z, fluid_cell_info_ptr);
}

void MpiMusic::GetHydroInfo_JETSCAPE(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  auto temp = bulk_info.get_tz(t, x, y, z);
  fluid_cell_info_ptr = std::unique_ptr<FluidCellInfo>(new FluidCellInfo(temp));
}

void MpiMusic::GetHydroInfo_MUSIC(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  fluid_cell_info_ptr = Jetscape::make_unique<FluidCellInfo>();
  fluidCell *fluidCell_ptr = new fluidCell;
  music_hydro_ptr->get_hydro_info(x, y, z, t, fluidCell_ptr);
  fluid_cell_info_ptr->energy_density = fluidCell_ptr->ed;
  fluid_cell_info_ptr->entropy_density = fluidCell_ptr->sd;
  fluid_cell_info_ptr->temperature = fluidCell_ptr->temperature;
  fluid_cell_info_ptr->pressure = fluidCell_ptr->pressure;
  fluid_cell_info_ptr->vx = fluidCell_ptr->vx;
  fluid_cell_info_ptr->vy = fluidCell_ptr->vy;
  fluid_cell_info_ptr->vz = fluidCell_ptr->vz;
  fluid_cell_info_ptr->mu_B = 0.0;
  fluid_cell_info_ptr->mu_C = 0.0;
  fluid_cell_info_ptr->mu_S = 0.0;
  fluid_cell_info_ptr->qgp_fraction = 0.0;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      fluid_cell_info_ptr->pi[i][j] = fluidCell_ptr->pi[i][j];
    }
  }
  fluid_cell_info_ptr->bulk_Pi = fluidCell_ptr->bulkPi;
  delete fluidCell_ptr;
}
