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
#include <vector>

#include "JetScapeLogger.h"
#include "CLViscWrapper.h"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<CLVisc> CLVisc::reg("CLVisc");

CLVisc::CLVisc() {
    hydro_status = NOT_START;
    SetId("CLVisc");
}


CLVisc::~CLVisc() {
}


void CLVisc::InitializeHydro(Parameter parameter_list) {
    JSINFO << "Initialize CLVisc ...";
    VERBOSE(8);

    std::string s = GetXMLElementText({"Hydro", "CLVisc", "name"});
    JSDEBUG << s << " to be initilizied ...";

    clvisc::Config cfg;

    std::string device_type = GetXMLElementText({"Hydro", "CLVisc", "device_type"});
        
    if(device_type == "cpu" || device_type == "CPU") {
        cfg.block_size = GetXMLElementInt({"Hydro", "CLVisc", "cpu_block_size"});
    } else {
        cfg.block_size = GetXMLElementInt({"Hydro", "CLVisc", "gpu_block_size"});
    }
    int device_id = GetXMLElementInt({"Hydro", "CLVisc", "device_id"});
    cfg.etaos_xmin = GetXMLElementInt({"Hydro", "CLVisc", "etaos_xmin"});
    cfg.etaos_ymin = GetXMLElementInt({"Hydro", "CLVisc", "etaos_ymin"});
    cfg.etaos_left_slop = GetXMLElementInt({"Hydro", "CLVisc", "etaos_left_slop"});
    cfg.etaos_right_slop = GetXMLElementInt({"Hydro", "CLVisc", "etaos_right_slop"});
    cfg.result_directory = GetXMLElementInt({"Hydro", "CLVisc", "result_directory"});
    doCooperFrye = GetXMLElementInt({"Hydro", "CLVisc", "Perform_CooperFrye_Feezeout"});
    cfg.tau0 = GetXMLElementDouble({"Hydro", "CLVisc", "tau0"});
    cfg.dt = GetXMLElementDouble({"Hydro", "CLVisc", "dtau"});
    cfg.ntskip = GetXMLElementDouble({"Hydro", "CLVisc", "ntau_skip"});
    cfg.nxskip = GetXMLElementDouble({"Hydro", "CLVisc", "nx_skip"});
    cfg.nyskip = GetXMLElementDouble({"Hydro", "CLVisc", "ny_skip"});
    cfg.nzskip = GetXMLElementDouble({"Hydro", "CLVisc", "netas_skip"});

    cfg.dx = ini->GetXStep();
    cfg.dy = ini->GetYStep();
    cfg.dz = ini->GetZStep();
    cfg.nx = ini->GetXSize();
    cfg.ny = ini->GetYSize();
    cfg.nz = ini->GetZSize();

    hydro_ = std::unique_ptr<clvisc::CLVisc>(new clvisc::CLVisc(cfg, device_type, device_id));
    initial_condition_scale_factor = GetXMLElementDouble({"Hydro", "CLVisc", "scale_factor"});
}

void CLVisc::EvolveHydro() {
    VERBOSE(8);
    JSINFO << "Initialize density profiles in CLVisc ...";
    std::vector<double> entropy_density = ini->GetEntropyDensityDistribution();
    double dx = ini->GetXStep();
    if (pre_eq_ptr == nullptr) {
        if (initial_condition_scale_factor == 1.0) {
            hydro_->read_ini(entropy_density);
        } else {
            std::for_each(entropy_density.begin(), entropy_density.end(),
                          [&](double & sd)
                          {sd = initial_condition_scale_factor * sd;});
            hydro_->read_ini(entropy_density);
        }
    } else {
        std::vector<double> vx_, vy_, vz_;
        for (size_t idx = 0; idx < pre_eq_ptr->ux_.size(); idx++) {
            vx_.push_back(pre_eq_ptr->ux_.at(idx)/pre_eq_ptr->utau_.at(idx));
            vy_.push_back(pre_eq_ptr->uy_.at(idx)/pre_eq_ptr->utau_.at(idx));
            vz_.push_back(pre_eq_ptr->ueta_.at(idx)/pre_eq_ptr->utau_.at(idx));
        }

        hydro_->read_ini(pre_eq_ptr->e_, vx_, vy_, vz_,
                         pre_eq_ptr->pi00_,
                         pre_eq_ptr->pi01_,
                         pre_eq_ptr->pi02_,
                         pre_eq_ptr->pi03_,
                         pre_eq_ptr->pi11_,
                         pre_eq_ptr->pi12_,
                         pre_eq_ptr->pi13_,
                         pre_eq_ptr->pi22_,
                         pre_eq_ptr->pi23_,
                         pre_eq_ptr->pi33_);
    }

    hydro_status = INITIALIZED;
    if (hydro_status == INITIALIZED) {
        JSINFO << "running CLVisc ...";
        hydro_->evolve();
        hydro_status = FINISHED;

        auto cfg = hydro_->get_config();
        float tau_min = cfg.tau0;
        float dtau = cfg.dt * cfg.ntskip;
        float dx = cfg.dx * cfg.nxskip;
        float dy = cfg.dy * cfg.nyskip;
        float detas = cfg.dz * cfg.nzskip;
        float x_min = - 0.5f * cfg.nx * cfg.dx;
        float y_min = - 0.5f * cfg.ny * cfg.dy;
        float etas_min = - 0.5f * cfg.nz * cfg.dz;
        int nx = int(floor((cfg.nx-1) / cfg.nxskip)) + 1;
        int ny = int(floor((cfg.ny-1) / cfg.nyskip)) + 1;
        int netas = int(floor((cfg.nz-1) / cfg.nzskip)) + 1;

        bulk_info.FromVector(hydro_->bulkinfo_.get_data(),
               hydro_->bulkinfo_.get_data_info(),
               tau_min, dtau, x_min, dx, nx,
               y_min, dy, ny, etas_min, detas, netas,
               false);
        //hydro_->bulkinfo_.save("bulk_data.csv");
    }
    if (hydro_status == FINISHED && doCooperFrye == 1) {
        JSINFO << "Cooper Frye not implemented yet";
    }
}

void CLVisc::GetHydroInfo(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
        std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr) {
      fluid_cell_info_ptr = make_unique<FluidCellInfo> ();
      if (hydro_status != FINISHED) {
        throw std::runtime_error("Hydro evolution is not finished ");
      }

      if (!bulk_info.tau_eta_is_tz) {
        Jetscape::real tau = std::sqrt(t * t - z * z);
        Jetscape::real eta = 0.5 * (std::log(t + z) - std::log(t - z));
        try {
            bulk_info.CheckInRange(tau, x, y, eta);
            *fluid_cell_info_ptr = bulk_info.get(tau, x, y, eta);
        } catch (std::exception& err) {
            JSWARN << err.what();
            //fluid_cell_info_ptr->Print();
        }
      } else {
        try {
            bulk_info.CheckInRange(t, x, y, z);
            *fluid_cell_info_ptr = bulk_info.get(t, x, y, z);
        } catch (std::exception & err) {
            JSWARN << err.what();
        }
      }
}
