// -----------------------------------------
// JETSCPAE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// Copyright [2018] <Chun Shen>
// -----------------------------------------

#include "SoftParticlization.h"

namespace Jetscape {

SoftParticlization::SoftParticlization() {
    Hadron_list_ = new std::vector< std::vector<Hadron>* >;
}

SoftParticlization::~SoftParticlization() {
    for (unsigned i = 0; i < Hadron_list_->size(); i++) {
        (*Hadron_list_)[i]->clear();
    }
    Hadron_list_->clear();
    delete Hadron_list_;
}

void SoftParticlization::Init() {
    JetScapeModuleBase::Init();
    INFO << "Intialize Soft particlization module ... " << GetId() << " ...";
    xml_ = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement(
                                                        "SoftParticlization");
    if (!xml_) {
        WARN << " : Missing XML SoftParticlization section in file!";
        exit(-1);
    }

    InitTask();
}

void SoftParticlization::Exec() {

}

void SoftParticlization::Clear() {

}

}  // end namespace Jetscape
