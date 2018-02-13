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
}

void iSS_CF::Clear() {
    VERBOSE(2) << "Finish the particle sampling";
    if (iSpectraSampler_ptr_ != nullptr) {
        delete iSpectraSampler_ptr_;
    }
}
