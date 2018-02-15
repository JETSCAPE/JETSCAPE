// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// Copyright [2018] <Chun Shen>
// -----------------------------------------

#include "JetScapeLogger.h"
#include "iSS_jetscape.h"

#include <string>

using namespace Jetscape;

iSS_CF::iSS_CF() {
    SetId("iSS");
    iSpectraSampler_ptr_ = nullptr;
}

iSS_CF::~iSS_CF() {
}

void iSS_CF::InitTask() {
    INFO << "Initialize a particle sampler (iSS)";
    iSS_xml_ = xml_->FirstChildElement("iSS");
    if (!iSS_xml_) {
        WARN << "No XML section for iSS! Please check the input file~";
        exit(-1);
    }
    string input_file = (
                iSS_xml_->FirstChildElement("iSS_input_file")->GetText());
    iSpectraSampler_ptr_ = new iSS;
    iSpectraSampler_ptr_->paraRdr_ptr->readFromFile(input_file);
    iSpectraSampler_ptr_->paraRdr_ptr->echo();
}

void iSS_CF::Exec() {
    int status = iSpectraSampler_ptr_->read_in_FO_surface();
    if (status != 0) {
        WARN << "Some errors happened in reading in the hyper-surface";
        exit(-1);
    }

    auto random_seed = (*get_mt19937_generator())();  // get random seed
    iSpectraSampler_ptr_->set_random_seed(random_seed);
    VERBOSE(2) << "Random seed used for the iSS module" << random_seed;
    
    status = iSpectraSampler_ptr_->generate_samples();
    if (status != 0) {
        WARN << "Some errors happened in generating particle samples";
        exit(-1);
    }
    pass_hadron_list_to_JETSCAPE();
}

void iSS_CF::Clear() {
    VERBOSE(2) << "Finish the particle sampling";
    if (iSpectraSampler_ptr_ != nullptr) {
        delete iSpectraSampler_ptr_;
    }
}

void iSS_CF::pass_hadron_list_to_JETSCAPE() {
    unsigned int nev = iSpectraSampler_ptr_->get_number_of_sampled_events();
    VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
    VERBOSE(4) << "number of events to pass : " << nev;
    for (unsigned int iev = 0; iev < nev; iev++) {
        Hadron_list_->push_back(new std::vector<Hadron>);
        unsigned int nparticles = (
                        iSpectraSampler_ptr_->get_number_of_particles(iev));
        VERBOSE(4) << "event " << iev << ": number of particles = "
                   << nparticles;
        for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
            iSS_Hadron current_hadron = (
                            iSpectraSampler_ptr_->get_hadron(iev, ipart));
            int hadron_label = 0;
            int hadron_status = -1;
            //int hadron_id = current_hadron.pid;
            int hadron_id = 1;   // just for testing need to be changed to the line above
            double hadron_mass = current_hadron.mass;
            FourVector hadron_p(current_hadron.px, current_hadron.py,
                                current_hadron.pz, current_hadron.E);
            FourVector hadron_x(current_hadron.x, current_hadron.y,
                                current_hadron.z, current_hadron.t);

            // create a JETSCAPE Hadron
            Hadron* jetscape_hadron = new Hadron(
                hadron_label, hadron_id, hadron_status, hadron_p, hadron_x);
            (*Hadron_list_)[iev]->push_back(*jetscape_hadron);
        }
    }
    VERBOSE(4) << "JETSCAPE received " << Hadron_list_->size() << " events.";
    for (unsigned int iev = 0; iev < Hadron_list_->size(); iev++) {
        VERBOSE(4) << "In event " << iev << " JETSCAPE received "
                   << (*Hadron_list_)[iev]->size() << " particles.";
    }
}
