// -----------------------------------------
// JETSCPAE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// Copyright [2018] <Chun Shen>
// -----------------------------------------

#include "SoftParticlization.h"

namespace Jetscape {

SoftParticlization::SoftParticlization() {

}

SoftParticlization::~SoftParticlization() {

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
}

void SoftParticlization::Exec() {

}

void SoftParticlization::Clear() {

}

}  // end namespace Jetscape