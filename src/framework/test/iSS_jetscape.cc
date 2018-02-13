// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// Copyright [2018] <Chun Shen>
// -----------------------------------------

#include "JetScapeLogger.h"
#include "iSS_jetscape.h"

using namespace Jetscape;

iSS_CF::iSS_CF() {
    SetId("iSS");
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
}

void iSS_CF::Exec() {
    auto random_seed = (*get_mt19937_generator())();
    VERBOSE(2) << "Random seed used for the iSS module" << random_seed;
}

void iSS_CF::Clear() {
    VERBOSE(2) << "Finish the particle sampling";
}
