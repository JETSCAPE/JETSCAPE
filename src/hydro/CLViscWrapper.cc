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

CLViscWrapper::CLViscWrapper() {
    hydro_status = NOT_START;
    doCooperFrye = 0;
    SetId("CLVisc");
}


CLViscWrapper::~CLViscWrapper() {
}


void CLViscWrapper::InitializeHydro(Parameter parameter_list) {
    INFO << "Initialize CLVisc ...";
    VERBOSE(8);
    tinyxml2::XMLElement *para =
                    GetHydroXML()->FirstChildElement("CLVisc");
    if (!para) {
        WARN << " : CLVisc not properly initialized in XML file ...";
        exit(-1);
    }

    clvisc::Config cfg;
    para->FirstChildElement("gpu_block_size")->QueryIntText(&cfg.block_size);
    std::string device_type = para->FirstChildElement("device_type")->GetText();
    std::cout << device_type << std::endl;
    int device_id;
    para->FirstChildElement("device_id")->QueryIntText(&device_id);
    para->FirstChildElement("etaos_xmin")->QueryFloatText(&cfg.etaos_xmin);
    para->FirstChildElement("etaos_ymin")->QueryFloatText(&cfg.etaos_ymin);
    para->FirstChildElement("etaos_left_slop")->QueryFloatText(&cfg.etaos_left_slop);
    para->FirstChildElement("etaos_right_slop")->QueryFloatText(&cfg.etaos_right_slop);
    cfg.result_directory = para->FirstChildElement("result_directory")->GetText();
    para->FirstChildElement("Perform_CooperFrye_Feezeout")->QueryIntText(&doCooperFrye);
    para->FirstChildElement("tau0")->QueryDoubleText(&cfg.tau0);
    para->FirstChildElement("dtau")->QueryDoubleText(&cfg.dt);

    cfg.dx = ini->GetXStep();
    cfg.dy = ini->GetYStep();
    cfg.dz = ini->GetZStep();
    cfg.nx = ini->GetXSize();
    cfg.ny = ini->GetYSize();
    cfg.nz = ini->GetZSize();

    hydro_ = std::unique_ptr<clvisc::CLVisc>(new clvisc::CLVisc(cfg, device_type, device_id));
    para->FirstChildElement("scale_factor")->QueryDoubleText(&initial_condition_scale_factor);
}

void CLViscWrapper::EvolveHydro() {
    VERBOSE(8);
    INFO << "Initialize density profiles in CLVisc ...";
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

        hydro_->read_ini(pre_eq_ptr->e_,
                         vx_,
                         vy_,
                         vz_,
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
        INFO << "running CLVisc ...";
        hydro_->evolve();
        hydro_status = FINISHED;
    }
    if (hydro_status == FINISHED && doCooperFrye == 1) {
        //music_hydro_ptr->run_Cooper_Frye(1);
        INFO << "Cooper Frye not implemented yet";
    }
}

void CLViscWrapper::GetHydroInfo(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
        std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr) {
    fluid_cell_info_ptr = std::make_unique<FluidCellInfo>();
}

