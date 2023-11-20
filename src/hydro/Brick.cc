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

#include <cstring>
#include <cmath>
#include <iostream>
#include <MakeUniqueHelper.h>

#include "JetScapeLogger.h"
#include "Brick.h"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<Brick> Brick::reg("Brick");

Brick::Brick() : FluidDynamics() {
  // initialize the parameter reader
  T_brick = 0.0; // GeV
  start_time = 0.0;
  bjorken_expansion_on = false;
  vx = 0;
  vy = 0;
  vz = 0;

  hydro_status = NOT_START;
  SetId("Brick");
  VERBOSE(8);
}

Brick::~Brick() { VERBOSE(8); }

void Brick::InitTask() {
  // kind of stupid ... do pointer GetHydroXML() via XML instance ...

  JSDEBUG << "Initialize Brick (Test) ...";
  VERBOSE(8);

  std::string s = GetXMLElementText({"Hydro", "Brick", "name"});
  JSDEBUG << s << " to be initialized ...";

  T_brick = GetXMLElementDouble({"Hydro", "Brick", "T"});
  JSDEBUG << s << " with T = " << T_brick;
  VERBOSE(2) << "Brick Temperature T = " << T_brick;

  tinyxml2::XMLElement *brick = GetXMLElement({"Hydro", "Brick"});
  if (brick->Attribute("bjorken_expansion_on", "true")) {
    bjorken_expansion_on = true;
    start_time = std::atof(brick->Attribute("start_time"));
  } else {
    if (brick->Attribute("start_time")){
      start_time = std::atof(brick->Attribute("start_time"));
    }
  }

  hydro_tau_0 = start_time;

  brick_L = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});

  //flow
  vx = GetXMLElementDouble({"Hydro", "Brick", "vx"});
  vy = GetXMLElementDouble({"Hydro", "Brick", "vy"});
  vz = GetXMLElementDouble({"Hydro", "Brick", "vz"});

  //Parameter parameter_list;
  GetParameterList().hydro_input_filename = (char *)"dummy"; //*(argv+1);
}

void Brick::InitializeHydro(Parameter parameter_list) {
  hydro_status = INITIALIZED;
}

void Brick::EvolveHydro() {
  VERBOSE(8);
  VERBOSE(2) << "size of sd = " << ini->GetEntropyDensityDistribution().size();
  hydro_status = FINISHED;
}

void Brick::GetHydroInfo(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    //                           FluidCellInfo* fluid_cell_info_ptr) {
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  // create the unique FluidCellInfo here
  fluid_cell_info_ptr = make_unique<FluidCellInfo>();

  // assign all the quantites to JETSCAPE output
  // thermodyanmic quantities

  if (hydro_status == FINISHED) {
    fluid_cell_info_ptr->energy_density = 0.0;
    fluid_cell_info_ptr->entropy_density = 0.0;
	if(t > brick_L){fluid_cell_info_ptr->temperature = 0.;}
    else if (bjorken_expansion_on) {
      fluid_cell_info_ptr->temperature =
          T_brick * std::pow(start_time / t, 1.0 / 3.0);
    } else {
      fluid_cell_info_ptr->temperature = T_brick;
    }
    fluid_cell_info_ptr->pressure = 0.0;
    // QGP fraction
    fluid_cell_info_ptr->qgp_fraction = 1.0;
    // chemical potentials
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    // dynamical quantites
    fluid_cell_info_ptr->vx = vx;
    fluid_cell_info_ptr->vy = vy;
    fluid_cell_info_ptr->vz = vz;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        fluid_cell_info_ptr->pi[i][j] = 0.0;
      }
    }
    fluid_cell_info_ptr->bulk_Pi = 0.0;
  } else {
    JSWARN << "Hydro not run yet ...";
    exit(-1);
  }
}
