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
#include "surfaceCell.h"

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
  int echoLevel = GetXMLElementInt({"vlevel"});
  music_hydro_ptr->set_parameter("JSechoLevel", echoLevel);

  flag_output_evo_to_file = (
      GetXMLElementInt({"Hydro", "MUSIC", "output_evolution_to_file"}));
  if (flag_output_evo_to_file == 1) {
    music_hydro_ptr->set_parameter("output_evolution_data", 2);
  } else {
    music_hydro_ptr->set_parameter("output_evolution_data", 0);
  }

  flag_output_evo_to_memory = (
      GetXMLElementInt({"Hydro", "MUSIC", "output_evolution_to_memory"}));
  if (flag_output_evo_to_memory == 1) {
    music_hydro_ptr->set_parameter("store_hydro_info_in_memory", 1);
    music_hydro_ptr->set_parameter("output_evolution_every_N_timesteps",
      GetXMLElementInt({"Hydro", "MUSIC", "output_evolution_every_N_timesteps"})
    );
  }

  double tau_hydro = (
          GetXMLElementDouble({"Hydro", "MUSIC", "Initial_time_tau_0"}));
  music_hydro_ptr->set_parameter("Initial_time_tau_0", tau_hydro);
  int beastMode = GetXMLElementInt({"Hydro", "MUSIC", "beastMode"});
  music_hydro_ptr->set_parameter("beastMode", beastMode);

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
    music_hydro_ptr->set_parameter("T_dependent_Shear_to_S_ratio",
                                   flag_shear_Tdep);
    if (flag_shear_Tdep == 3) {
      double shear_kinkT = (
        GetXMLElementDouble({"Hydro", "MUSIC", "eta_over_s_T_kink_in_GeV"}));
      music_hydro_ptr->set_parameter("shear_viscosity_3_T_kink_in_GeV",
                                     shear_kinkT);
      double shear_lowTslope = (
        GetXMLElementDouble({"Hydro", "MUSIC",
                             "eta_over_s_low_T_slope_in_GeV"}));
      music_hydro_ptr->set_parameter("shear_viscosity_3_low_T_slope_in_GeV",
                                     shear_lowTslope);
      double shear_highTslope = (
        GetXMLElementDouble({"Hydro", "MUSIC",
                             "eta_over_s_high_T_slope_in_GeV"}));
      music_hydro_ptr->set_parameter("shear_viscosity_3_high_T_slope_in_GeV",
                                     shear_highTslope);
      double shear_kink = (
        GetXMLElementDouble({"Hydro", "MUSIC", "eta_over_s_at_kink"}));
      music_hydro_ptr->set_parameter("shear_viscosity_3_at_kink", shear_kink);
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
        music_hydro_ptr->set_parameter("bulk_viscosity_3_max", bulk_max);
        double bulk_peakT = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_T_peak_in_GeV"});
        music_hydro_ptr->set_parameter("bulk_viscosity_3_T_peak_in_GeV",
                                       bulk_peakT);
        double bulk_width = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_width_in_GeV"});
        music_hydro_ptr->set_parameter("bulk_viscosity_3_width_in_GeV",
                                       bulk_width);
        double bulk_asy = GetXMLElementDouble(
              {"Hydro", "MUSIC", "zeta_over_s_lambda_asymm"});
        music_hydro_ptr->set_parameter("bulk_viscosity_3_lambda_asymm",
                                       bulk_asy);
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

  flag_surface_in_memory = GetXMLElementInt(
          {"Hydro", "MUSIC", "surface_in_memory"});
  if (flag_surface_in_memory == 1) {
    music_hydro_ptr->set_parameter("surface_in_memory", 1);
  } else {
    music_hydro_ptr->set_parameter("surface_in_memory", 0);
  }

  music_hydro_ptr->add_hydro_source_terms(hydro_source_terms_ptr);
}

void MpiMusic::EvolveHydro() {
  VERBOSE(8);
  JSINFO << "Initialize density profiles in MUSIC ...";

  if (pre_eq_ptr == nullptr) {
    JSWARN << "Missing the pre-equilibrium module ...";
    exit(1);
  } else {
    double dx = ini->GetXStep();
    double dz = ini->GetZStep();
    double z_max = ini->GetZMax();
    int nz = ini->GetZSize();
    double tau0 = pre_eq_ptr->GetPreequilibriumEndTime();
    music_hydro_ptr->initialize_hydro_from_jetscape_preequilibrium_vectors(
        tau0,
        dx, dz, z_max, nz, pre_eq_ptr->e_, pre_eq_ptr->P_,
        pre_eq_ptr->utau_, pre_eq_ptr->ux_, pre_eq_ptr->uy_, pre_eq_ptr->ueta_,
        pre_eq_ptr->pi00_, pre_eq_ptr->pi01_, pre_eq_ptr->pi02_,
        pre_eq_ptr->pi03_, pre_eq_ptr->pi11_, pre_eq_ptr->pi12_,
        pre_eq_ptr->pi13_, pre_eq_ptr->pi22_, pre_eq_ptr->pi23_,
        pre_eq_ptr->pi33_, pre_eq_ptr->bulk_Pi_);
    JSINFO << "initial density profile dx = " << dx << " fm";
  }

  has_source_terms = false;
  if (hydro_source_terms_ptr->get_number_of_sources() > 0) {
    has_source_terms = true;
  }
  JSINFO << "number of source terms: "
         << hydro_source_terms_ptr->get_number_of_sources()
         << ", total E = " << hydro_source_terms_ptr->get_total_E_of_sources()
         << " GeV.";

  if (pre_eq_ptr != nullptr) {
    double tau0 = pre_eq_ptr->GetPreequilibriumEndTime();
    JSINFO << "hydro initial time set by PreEq module tau0 = "
           << tau0 << " fm/c";
    if (flag_output_evo_to_memory == 1) {
      // need to ensure preEq and hydro use the same dtau so that
      // the combined evolution history file is properly set
      double dtau = pre_eq_ptr->GetPreequilibriumEvodtau();
      JSINFO << "Reset MUSIC dtau by PreEq module: dtau = " << dtau << " fm/c";
      music_hydro_ptr->set_parameter("dtau", dtau);
      if (!has_source_terms) {
        // only the preEq evo with the first hydro without source term
        // will be stored in memory for jet energy loss calculations
        clear_up_evolution_data();
        PassPreEqEvolutionHistoryToFramework();
      }
    }
  }

  hydro_status = INITIALIZED;

  if (hydro_status == INITIALIZED) {
    JSINFO << "running MUSIC ...";
    music_hydro_ptr->run_hydro();
    hydro_status = FINISHED;
  }

  if (flag_output_evo_to_memory == 1) {
    if (!has_source_terms) {
      // only the first hydro without source term will be stored
      // in memory for jet energy loss calculations
      if (pre_eq_ptr == nullptr) {
        clear_up_evolution_data();
      }
      PassHydroEvolutionHistoryToFramework();
      JSINFO << "number of fluid cells received by the JETSCAPE: "
             << bulk_info.data.size();
    }
    music_hydro_ptr->clear_hydro_info_from_memory();
  }

  if (flag_output_evo_to_file == 1) {
    //if (!has_source_terms) {
    //  // only the first hydro without source term will be stored
    //  // in memory for jet energy loss calculations
    //  PassHydroEvolutionHistoryToFramework();
    //  JSINFO << "number of fluid cells received by the JETSCAPE: "
    //         << bulk_info.data.size();
    //}
    //music_hydro_ptr->clear_hydro_info_from_memory();

    // add hydro_id to the hydro evolution filename
    std::ostringstream system_command;
    system_command << "mv evolution_all_xyeta.dat "
                   << "evolution_all_xyeta_" << GetId() << ".dat";
    system(system_command.str().c_str());

    //std::vector<SurfaceCellInfo> surface_cells;
    //if (freezeout_temperature > 0.0) {
    //  FindAConstantTemperatureSurface(freezeout_temperature, surface_cells);
    //}
  }

  if (flag_surface_in_memory == 1) {
    clearSurfaceCellVector();
    PassHydroSurfaceToFramework();
  } else {
    collect_freeze_out_surface();
  }

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


void MpiMusic::SetPreEqGridInfo() {
  bulk_info.tau_min = pre_eq_ptr->GetPreequilibriumStartTime();
  bulk_info.dtau = pre_eq_ptr->GetPreequilibriumEvodtau();
  JSINFO << "preEq evo: tau_0 = " << bulk_info.tau_min
         << " fm/c, dtau = " << bulk_info.dtau << " fm/c.";
}


void MpiMusic::SetHydroGridInfo() {
  bulk_info.neta = music_hydro_ptr->get_neta();
  bulk_info.nx = music_hydro_ptr->get_nx();
  bulk_info.ny = music_hydro_ptr->get_nx();
  bulk_info.x_min = -music_hydro_ptr->get_hydro_x_max();
  bulk_info.dx = music_hydro_ptr->get_hydro_dx();
  bulk_info.y_min = -music_hydro_ptr->get_hydro_x_max();
  bulk_info.dy = music_hydro_ptr->get_hydro_dx();
  bulk_info.eta_min = -music_hydro_ptr->get_hydro_eta_max();
  bulk_info.deta = music_hydro_ptr->get_hydro_deta();

  bulk_info.boost_invariant = music_hydro_ptr->is_boost_invariant();

  if (pre_eq_ptr == nullptr) {
    bulk_info.tau_min = music_hydro_ptr->get_hydro_tau0();
    bulk_info.dtau = music_hydro_ptr->get_hydro_dtau();
    bulk_info.ntau = music_hydro_ptr->get_ntau();
  } else {
    bulk_info.ntau = music_hydro_ptr->get_ntau() + pre_eq_ptr->get_ntau();
  }
}


void MpiMusic::PassPreEqEvolutionHistoryToFramework() {
  JSINFO << "Passing preEq evolution information to JETSCAPE ... ";
  auto number_of_cells = pre_eq_ptr->get_number_of_fluid_cells();
  JSINFO << "total number of preEq fluid cells: " << number_of_cells;

  SetPreEqGridInfo();

  for (int i = 0; i < number_of_cells; i++) {
    std::unique_ptr<FluidCellInfo> fluid_cell_info_ptr(new FluidCellInfo);
    pre_eq_ptr->get_fluid_cell_with_index(i, fluid_cell_info_ptr);
    StoreHydroEvolutionHistory(fluid_cell_info_ptr);
  }
  pre_eq_ptr->clear_evolution_data();
}


void MpiMusic::PassHydroEvolutionHistoryToFramework() {
  JSINFO << "Passing hydro evolution information to JETSCAPE ... ";
  auto number_of_cells = music_hydro_ptr->get_number_of_fluid_cells();
  JSINFO << "total number of MUSIC fluid cells: " << number_of_cells;

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


void MpiMusic::PassHydroSurfaceToFramework() {
    JSINFO << "Passing hydro surface cells to JETSCAPE ... ";
    auto number_of_cells = music_hydro_ptr->get_number_of_surface_cells();
    JSINFO << "total number of fluid cells: " << number_of_cells;
    SurfaceCell surfaceCell_i;
    for (int i = 0; i < number_of_cells; i++) {
        SurfaceCellInfo surface_cell_info;
        music_hydro_ptr->get_surface_cell_with_index(i, surfaceCell_i);
        surface_cell_info.tau = surfaceCell_i.xmu[0];
        surface_cell_info.x = surfaceCell_i.xmu[1];
        surface_cell_info.y = surfaceCell_i.xmu[2];
        surface_cell_info.eta = surfaceCell_i.xmu[3];
        double u[4];
        for (int j = 0; j < 4; j++) {
            surface_cell_info.d3sigma_mu[j] = surfaceCell_i.d3sigma_mu[j];
            surface_cell_info.umu[j] = surfaceCell_i.umu[j];
        }
        surface_cell_info.energy_density = surfaceCell_i.energy_density;
        surface_cell_info.temperature = surfaceCell_i.temperature;
        surface_cell_info.pressure = surfaceCell_i.pressure;
        surface_cell_info.baryon_density = surfaceCell_i.rho_b;
        surface_cell_info.mu_B = surfaceCell_i.mu_B;
        surface_cell_info.mu_Q = surfaceCell_i.mu_Q;
        surface_cell_info.mu_S = surfaceCell_i.mu_S;
        for (int j = 0; j < 10; j++) {
            surface_cell_info.pi[j] = surfaceCell_i.shear_pi[j];
        }
        surface_cell_info.bulk_Pi = surfaceCell_i.bulk_Pi;
        StoreSurfaceCell(surface_cell_info);
    }
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
