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
// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// -----------------------------------------

#include "JetScapeLogger.h"
#include "iSpectraSamplerWrapper.h"

#include <memory>
#include <string>
#include <fstream>

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<iSpectraSamplerWrapper>
    iSpectraSamplerWrapper::reg("iSS");

iSpectraSamplerWrapper::iSpectraSamplerWrapper() {
    SetId("iSS");
    statusCode_ = 0;
}

iSpectraSamplerWrapper::~iSpectraSamplerWrapper() {}

void iSpectraSamplerWrapper::InitTask() {

  JSINFO << "Initialize a particle sampler (iSS)";

  std::string input_file =
      GetXMLElementText({"SoftParticlization", "iSS", "iSS_input_file"});
  std::string table_path =
      GetXMLElementText({"SoftParticlization", "iSS", "iSS_table_path"});
  std::string particle_table_path =
      GetXMLElementText({"SoftParticlization", "iSS",
                         "iSS_particle_table_path"});
  std::string working_path =
      GetXMLElementText({"SoftParticlization", "iSS", "iSS_working_path"});
  int hydro_mode =
      GetXMLElementInt({"SoftParticlization", "iSS", "hydro_mode"});
  int number_of_repeated_sampling = GetXMLElementInt(
      {"SoftParticlization", "iSS", "number_of_repeated_sampling"});
  int flag_perform_decays = GetXMLElementInt(
      {"SoftParticlization", "iSS", "Perform_resonance_decays"});
  int afterburner_type = (
      GetXMLElementInt({"SoftParticlization", "iSS", "afterburner_type"}));

  if (!boost_invariance) {
    hydro_mode = 2;
  }

  iSpectraSampler_ptr_ = std::unique_ptr<iSS>(
          new iSS(working_path, table_path, particle_table_path, input_file));
  iSpectraSampler_ptr_->paraRdr_ptr->readFromFile(input_file);

  // overwrite some parameters
  int echoLevel = GetXMLElementInt({"vlevel"});
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("JSechoLevel", echoLevel);

  iSpectraSampler_ptr_->paraRdr_ptr->setVal("hydro_mode", hydro_mode);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("afterburner_type",
                                            afterburner_type);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("output_samples_into_files", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_OSCAR_format", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_gzip_format", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_binary_format", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("store_samples_in_memory", 1);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("number_of_repeated_sampling",
                                            number_of_repeated_sampling);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("perform_decays",
                                            flag_perform_decays);

  // set default parameters
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("turn_on_shear", 1);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("turn_on_bulk", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("turn_on_rhob", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("turn_on_diff", 0);


  iSpectraSampler_ptr_->paraRdr_ptr->setVal("restrict_deltaf", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("deltaf_max_ratio", 1.0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("f0_is_not_small", 1);

  iSpectraSampler_ptr_->paraRdr_ptr->setVal("calculate_vn", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("MC_sampling", 4);
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_spectators", 0);
 
  iSpectraSampler_ptr_->paraRdr_ptr->setVal("RegVisYield", 1);

  iSpectraSampler_ptr_->paraRdr_ptr->setVal(
      "sample_upto_desired_particle_number", 0);
  iSpectraSampler_ptr_->paraRdr_ptr->echo();
}

void iSpectraSamplerWrapper::Exec() {
  JSINFO << "running iSS ...";

  // generate symbolic links with music_input_file
  std::string music_input_file_path = GetXMLElementText(
          {"Hydro", "MUSIC", "MUSIC_input_file"});
  std::string working_path =
      GetXMLElementText({"SoftParticlization", "iSS", "iSS_working_path"});
  std::string music_input = working_path + "/music_input";
  std::ifstream inputfile(music_input.c_str());
  if (!inputfile.good()) {
    std::ostringstream system_command;
    system_command << "ln -s " << music_input_file_path << " "
                   << music_input;
    system(system_command.str().c_str());
  }
  inputfile.close();

  int nCells = getSurfCellVector();
  if (nCells == 0) {
    int status = iSpectraSampler_ptr_->read_in_FO_surface();
    if (status != 0) {
      JSWARN << "Some errors happened in reading in the hyper-surface";
      exit(-1);
    }
    nCells = 1;
  }

  auto random_seed = (*GetMt19937Generator())(); // get random seed
  iSpectraSampler_ptr_->set_random_seed(random_seed);
  VERBOSE(2) << "Random seed used for the iSS module" << random_seed;

  if (nCells > 0) {
    int status = iSpectraSampler_ptr_->generate_samples();
    if (status != 0) {
      JSWARN << "Some errors happened in generating particle samples";
      exit(-1);
    }
    PassHadronListToJetscape();
  }
  JSINFO << "iSS finished.";
}

void iSpectraSamplerWrapper::Clear() {
  VERBOSE(2) << "Finish the particle sampling";
  iSpectraSampler_ptr_->clear();
  for (unsigned i = 0; i < Hadron_list_.size(); i++) {
    Hadron_list_.at(i).clear();
  }
  Hadron_list_.clear();
}

void iSpectraSamplerWrapper::PassHadronListToJetscape() {
  unsigned int nev = iSpectraSampler_ptr_->get_number_of_sampled_events();
  VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
  VERBOSE(4) << "number of events to pass : " << nev;
  for (unsigned int iev = 0; iev < nev; iev++) {
    std::vector<shared_ptr<Hadron>> hadrons;
    unsigned int nparticles =
        (iSpectraSampler_ptr_->get_number_of_particles(iev));
    VERBOSE(4) << "event " << iev << ": number of particles = " << nparticles;
    for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
      iSS_Hadron current_hadron =
          (iSpectraSampler_ptr_->get_hadron(iev, ipart));
      int hadron_label = 0;
      int hadron_status = 11;
      int hadron_id = current_hadron.pid;
      //int hadron_id = 1;   // just for testing need to be changed to the line above
      double hadron_mass = current_hadron.mass;
      FourVector hadron_p(current_hadron.px, current_hadron.py,
                          current_hadron.pz, current_hadron.E);
      FourVector hadron_x(current_hadron.x, current_hadron.y, current_hadron.z,
                          current_hadron.t);

      // create a JETSCAPE Hadron
      hadrons.push_back(make_shared<Hadron>(hadron_label, hadron_id,
                                            hadron_status, hadron_p, hadron_x,
                                            hadron_mass));
      //Hadron* jetscape_hadron = new Hadron(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass);
      //(*Hadron_list_)[iev]->push_back(*jetscape_hadron);
    }
    Hadron_list_.push_back(hadrons);
  }
  if (nev > 0) {
    VERBOSE(4) << "JETSCAPE received " << Hadron_list_.size() << " events.";
    for (unsigned int iev = 0; iev < Hadron_list_.size(); iev++) {
      VERBOSE(4) << "In event " << iev << " JETSCAPE received "
                 << Hadron_list_.at(iev).size() << " particles.";
    }
  }
}

void iSpectraSamplerWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(4) << "In iSpectraSamplerWrapper::WriteTask";
  auto f = w.lock();
  if (!f)
    return;

  f->WriteComment("JetScape module: " + GetId());
  if (Hadron_list_.size() > 0) {
    f->WriteComment("Final State Bulk Hadrons");
    for (unsigned int j = 0; j < Hadron_list_.size(); j++) {
      vector<shared_ptr<Hadron>> hadVec = Hadron_list_.at(j);
      for (unsigned int i = 0; i < hadVec.size(); i++) {
        f->WriteWhiteSpace("[" + to_string(i) + "] H");
        f->Write(hadVec.at(i));
      }
    }
  } else {
    f->WriteComment("There are no bulk Hadrons");
  }
}


int iSpectraSamplerWrapper::getSurfCellVector() {
  std::vector<SurfaceCellInfo> surfVec;
  std::vector<FO_surf> FOsurf_array;
  GetHydroHyperSurface(surfVec);
  int nCells = surfVec.size();
  JSINFO << "surface cell size: " << nCells;
  for (const auto surf_i: surfVec) {
    FO_surf iSS_surf_cell;
    iSS_surf_cell.tau = surf_i.tau;
    iSS_surf_cell.xpt = surf_i.x;
    iSS_surf_cell.ypt = surf_i.y;
    iSS_surf_cell.eta = surf_i.eta;
    iSS_surf_cell.da0 = surf_i.d3sigma_mu[0];
    iSS_surf_cell.da1 = surf_i.d3sigma_mu[1];
    iSS_surf_cell.da2 = surf_i.d3sigma_mu[2];
    iSS_surf_cell.da3 = surf_i.d3sigma_mu[3];
    iSS_surf_cell.u0 = surf_i.umu[0];
    iSS_surf_cell.u1 = surf_i.umu[1];
    iSS_surf_cell.u2 = surf_i.umu[2];
    iSS_surf_cell.u3 = surf_i.umu[3];
    iSS_surf_cell.Edec = surf_i.energy_density;
    iSS_surf_cell.Tdec = surf_i.temperature;
    iSS_surf_cell.Pdec = surf_i.pressure;
    iSS_surf_cell.muB = surf_i.mu_B;
    iSS_surf_cell.muQ = surf_i.mu_Q;
    iSS_surf_cell.muS = surf_i.mu_S;
    iSS_surf_cell.pi00 = surf_i.pi[0];
    iSS_surf_cell.pi01 = surf_i.pi[1];
    iSS_surf_cell.pi02 = surf_i.pi[2];
    iSS_surf_cell.pi03 = surf_i.pi[3];
    iSS_surf_cell.pi11 = surf_i.pi[4];
    iSS_surf_cell.pi12 = surf_i.pi[5];
    iSS_surf_cell.pi13 = surf_i.pi[6];
    iSS_surf_cell.pi22 = surf_i.pi[7];
    iSS_surf_cell.pi23 = surf_i.pi[8];
    iSS_surf_cell.pi33 = surf_i.pi[9];
    iSS_surf_cell.bulkPi = surf_i.bulk_Pi;
    FOsurf_array.push_back(iSS_surf_cell);
  }
  iSpectraSampler_ptr_->getSurfaceCellFromJETSCAPE(FOsurf_array);
  ClearHydroHyperSurface();
  return(nCells);
}
