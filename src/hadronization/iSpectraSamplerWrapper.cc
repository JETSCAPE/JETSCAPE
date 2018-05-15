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

#include <string>

using namespace Jetscape;

iSpectraSamplerWrapper::iSpectraSamplerWrapper() {
    SetId("iSS");
    iSpectraSampler_ptr_ = nullptr;
}

iSpectraSamplerWrapper::~iSpectraSamplerWrapper() {
}

void iSpectraSamplerWrapper::InitTask() {
    INFO << "Initialize a particle sampler (iSS)";
    iSS_xml_ = xml_->FirstChildElement("iSS");
    if (!iSS_xml_) {
        WARN << "No XML section for iSS! Please check the input file~";
        exit(-1);
    }
    string input_file = (
                iSS_xml_->FirstChildElement("iSS_input_file")->GetText());
    string working_path = (
                iSS_xml_->FirstChildElement("iSS_working_path")->GetText());
    int hydro_mode;
    iSS_xml_->FirstChildElement("hydro_mode")->QueryIntText(&hydro_mode);

    int number_of_repeated_sampling;
    iSS_xml_->FirstChildElement("number_of_repeated_sampling")->QueryIntText(
                                            &number_of_repeated_sampling);

    int flag_perform_decays;
    iSS_xml_->FirstChildElement("Perform_resonance_decays")->QueryIntText(
                                            &flag_perform_decays);

    iSpectraSampler_ptr_ = new iSS(working_path);
    iSpectraSampler_ptr_->paraRdr_ptr->readFromFile(input_file);

    // overwrite some parameters
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("hydro_mode", hydro_mode);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("output_samples_into_files", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_OSCAR_format", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("use_gzip_format", 0);
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
    
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_shear", 1);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_bulk", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("bulk_deltaf_kind", 1);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("include_deltaf_diffusion", 0);
    
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("restrict_deltaf", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("deltaf_max_ratio", 1.0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("f0_is_not_small", 1);

    iSpectraSampler_ptr_->paraRdr_ptr->setVal("calculate_vn", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->setVal("MC_sampling", 2);

    iSpectraSampler_ptr_->paraRdr_ptr->setVal(
                                    "sample_upto_desired_particle_number", 0);
    iSpectraSampler_ptr_->paraRdr_ptr->echo();
}

void iSpectraSamplerWrapper::Exec() {
    int status = iSpectraSampler_ptr_->read_in_FO_surface();
    if (status != 0) {
        WARN << "Some errors happened in reading in the hyper-surface";
        exit(-1);
    }

    auto random_seed = (*GetMt19937Generator())();  // get random seed
    iSpectraSampler_ptr_->set_random_seed(random_seed);
    VERBOSE(2) << "Random seed used for the iSS module" << random_seed;
    
    status = iSpectraSampler_ptr_->generate_samples();
    if (status != 0) {
        WARN << "Some errors happened in generating particle samples";
        exit(-1);
    }
    PassHadronListToJetscape();
}

void iSpectraSamplerWrapper::Clear() {
    VERBOSE(2) << "Finish the particle sampling";
    if (iSpectraSampler_ptr_ != nullptr) {
        delete iSpectraSampler_ptr_;
    }
}

void iSpectraSamplerWrapper::PassHadronListToJetscape() {
    unsigned int nev = iSpectraSampler_ptr_->get_number_of_sampled_events();
    VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
    VERBOSE(4) << "number of events to pass : " << nev;
    for (unsigned int iev = 0; iev < nev; iev++) {
        std::vector<shared_ptr<Hadron>> hadrons;
        unsigned int nparticles = (iSpectraSampler_ptr_->get_number_of_particles(iev));
        VERBOSE(4) << "event " << iev << ": number of particles = "<< nparticles;
        for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
            iSS_Hadron current_hadron = (iSpectraSampler_ptr_->get_hadron(iev, ipart));
            int hadron_label = 0;
            int hadron_status = -1;
            int hadron_id = current_hadron.pid;
            //int hadron_id = 1;   // just for testing need to be changed to the line above
            double hadron_mass = current_hadron.mass;
            FourVector hadron_p(current_hadron.px, current_hadron.py,
                                current_hadron.pz, current_hadron.E);
            FourVector hadron_x(current_hadron.x, current_hadron.y,
                                current_hadron.z, current_hadron.t);

            // create a JETSCAPE Hadron
            hadrons.push_back(make_shared<Hadron>(hadron_label,hadron_id,hadron_status,hadron_p,hadron_x, hadron_mass));
            //Hadron* jetscape_hadron = new Hadron(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass);
            //(*Hadron_list_)[iev]->push_back(*jetscape_hadron);
        }
        Hadron_list_.push_back(hadrons);
    }
    VERBOSE(4) << "JETSCAPE received " << Hadron_list_.size() << " events.";
    for (unsigned int iev = 0; iev < Hadron_list_.size(); iev++) {
        VERBOSE(4) << "In event " << iev << " JETSCAPE received " << Hadron_list_.at(iev).size() << " particles.";
    }

}

void iSpectraSamplerWrapper::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(4)<<"In iSpectraSamplerWrapper::WriteTask";
  auto f = w.lock();
  if ( !f ) return;

  f->WriteComment("JetScape module: "+GetId());
  if(Hadron_list_.size()>0) {
    f->WriteComment("Final State Bulk Hadrons");
    for(unsigned int j=0; j<Hadron_list_.size(); j++){
      vector<shared_ptr<Hadron>> hadVec = Hadron_list_.at(j);
      for(unsigned int i=0; i<hadVec.size(); i++) {
	f->WriteWhiteSpace("["+to_string(i)+"] H");
	f->Write(hadVec.at(i));
      }
    }
  } else {
    f->WriteComment("There are no bulk Hadrons");
  }
  
}

