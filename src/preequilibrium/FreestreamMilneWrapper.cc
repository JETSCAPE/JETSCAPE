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

#include "JetScapeLogger.h"
#include "FreestreamMilneWrapper.h"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<FreestreamMilneWrapper>
    FreestreamMilneWrapper::reg("FreestreamMilne");

FreestreamMilneWrapper::FreestreamMilneWrapper() {
  preequilibrium_status_ = NOT_STARTED;
  SetId("Freestream-Milne");
}

FreestreamMilneWrapper::~FreestreamMilneWrapper() {
  if (preequilibrium_status_ != NOT_STARTED)
    delete fsmilne_ptr;
}

void FreestreamMilneWrapper::InitializePreequilibrium(
    PreEquilibriumParameterFile parameter_list) {
  JSINFO << "Initialize freestream-milne ...";
  VERBOSE(8);

  std::string input_file = GetXMLElementText(
      {"Preequilibrium", "FreestreamMilne", "freestream_input_file"});
  //is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory

  fsmilne_ptr = new FREESTREAMMILNE();
  struct parameters *params = fsmilne_ptr->configure(input_file.c_str());

  //overwriting tau0,tauj,taus from xml file
  //tau0: initial time for initial condition
  //tauj: initial output time of background for hard probe
  //taus: end time for freestream or initial time for hydro
  //dtau: the free-streaming time
  double tau0 = GetXMLElementDouble(
      {"Preequilibrium", "tau0"});
  double tauj = GetXMLElementDouble(
      {"Preequilibrium", "tauj"});
  double taus = GetXMLElementDouble(
      {"Preequilibrium", "taus"});
  int FlagEvo = GetXMLElementDouble(
      {"Preequilibrium", "evolutionInMemory"});

  params->TAU0 = tau0;
  params->TAUJ = tauj;
  params->DTAU = taus - tau0;
  params->evolutionInMemory = FlagEvo;

  //settings for the grid size
  int nx = ini->GetXSize();
  int ny = ini->GetYSize();
  int neta = ini->GetZSize();
  params->DIM_X = nx;
  params->DIM_Y = ny;
  params->DIM_ETA = neta;

  //settings for the grid step size
  double dx = ini->GetXStep();
  double dy = ini->GetYStep();
  double deta = ini->GetZStep();
  params->DX = dx;
  params->DY = dy;
  params->DETA = deta;

  //setting for the number of time steps
  int ntau = GetXMLElementInt({"Preequilibrium", "FreestreamMilne", "ntau"});
  params->NT = ntau;
}

void FreestreamMilneWrapper::EvolvePreequilibrium() {
  VERBOSE(8);
  JSINFO << "Initialize energy density profile in freestream-milne ...";
  // grab initial energy density from vector from initial state module
  std::vector<double> entropy_density =
      ini->GetEntropyDensityDistribution(); //note that this is the energy density when read by freestream-milne, not actually the entropy density!
  std::vector<float> entropy_density_float(entropy_density.begin(),
                                           entropy_density.end());
  fsmilne_ptr->initialize_from_vector(entropy_density_float);
  JSINFO << " TRENTO event generated and loaded ";


  preequilibrium_status_ = INIT;
  if (preequilibrium_status_ == INIT) {
    JSINFO << "running freestream-milne ...";
    // evolve the medium via freestreaming
    fsmilne_ptr->run_freestream_milne();
    preequilibrium_status_ = DONE;
  }

  // now prepare to send the resulting hydro variables to the hydro module by coping hydro vectors to Preequilibrium base class members
 preequilibrium_tau_max_ = fsmilne_ptr->tau_LandauMatch;
  fsmilne_ptr->output_to_vectors(e_, P_, utau_, ux_, uy_, ueta_, pi00_, pi01_,
                                 pi02_, pi03_, pi11_, pi12_, pi13_, pi22_,
                                 pi23_, pi33_, bulk_Pi_);
}


void FreestreamMilneWrapper::get_fluid_cell_with_index(
        const int idx, std::unique_ptr<FluidCellInfo> &info_ptr) {
    fluidCell fluidCell_ptr;
    fsmilne_ptr->get_fluid_cell_with_index(idx, fluidCell_ptr);
    info_ptr->energy_density = fluidCell_ptr.ed;
    info_ptr->entropy_density = fluidCell_ptr.sd;
    info_ptr->temperature = fluidCell_ptr.temperature;
    info_ptr->pressure = fluidCell_ptr.pressure;
    info_ptr->vx = fluidCell_ptr.vx;
    info_ptr->vy = fluidCell_ptr.vy;
    info_ptr->vz = fluidCell_ptr.vz;
    info_ptr->mu_B = 0.0;
    info_ptr->mu_C = 0.0;
    info_ptr->mu_S = 0.0;
    info_ptr->qgp_fraction = 0.0;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        info_ptr->pi[i][j] = fluidCell_ptr.pi[i][j];
      }
    }
    info_ptr->bulk_Pi = fluidCell_ptr.bulkPi;
}
