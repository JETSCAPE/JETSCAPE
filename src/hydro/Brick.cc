/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
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

Brick::Brick() : FluidDynamics() {
    // initialize the parameter reader
    T_brick = 0.0;  // GeV
    start_time = 0.0;
    bjorken_expansion_on = false;

    hydro_status = NOT_START;
    SetId("Brick");
    VERBOSE(8);
}


Brick::~Brick() {VERBOSE(8);}

void Brick::InitTask()
{
  // kind of stupid ... do pointer GetHydroXML() via XML instance ...
  
  JSDEBUG<<"Initialize Brick (Test) ...";
  VERBOSE(8);
  tinyxml2::XMLElement *brick=GetHydroXML()->FirstChildElement("Brick");

  if (brick) {
      string s = brick->FirstChildElement( "name" )->GetText();

      JSDEBUG << s << " to be initilizied ...";
      
      brick->FirstChildElement("T")->QueryDoubleText(&T_brick);

      JSDEBUG << s << " with T = "<<T_brick;
      VERBOSE(2)<<"Brick Temperature T = "<<T_brick;

      if ( brick->Attribute("bjorken_expansion_on", "true") ) {
          bjorken_expansion_on = true;
          start_time = std::atof(brick->Attribute("start_time"));
      }

      //Parameter parameter_list;
      GetParameterList().hydro_input_filename = (char*) "dummy"; //*(argv+1);

    } else {
      WARN << " : Brick not properly initialized in XML file ...";
      exit(-1);
    }
}

void Brick::InitializeHydro(Parameter parameter_list) {
  hydro_status = INITIALIZED;
}


void Brick::EvolveHydro() {
  VERBOSE(8);
  VERBOSE(2) << "size of sd = " << ini->GetEntropyDensityDistribution().size();
  hydro_status = FINISHED;
}


void Brick::GetHydroInfo(real t, real x, real y, real z,
			   //                           FluidCellInfo* fluid_cell_info_ptr) {
			   std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr){
    // create the unique FluidCellInfo here
    fluid_cell_info_ptr=std::make_unique<FluidCellInfo>();

    // assign all the quantites to JETSCAPE output
    // thermodyanmic quantities

  if (hydro_status == FINISHED)
    {
      fluid_cell_info_ptr->energy_density = 0.0;
      fluid_cell_info_ptr->entropy_density = 0.0;
      if ( bjorken_expansion_on ) {
          fluid_cell_info_ptr->temperature = T_brick * std::pow(start_time/t, 1.0/3.0);
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
  else
    {
      WARN<<"Hydro not run yet ...";
      exit(-1);
    }
}
