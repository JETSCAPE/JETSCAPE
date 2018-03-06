// -----------------------------------------
// JETSCPAE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// Copyright [2018] <Chun Shen>
// -----------------------------------------


#ifndef SOFTPARTICLIZATION_H_
#define SOFTPARTICLIZATION_H_

#include <vector>

#include "JetScapeModuleBase.h"
#include "JetClass.hpp"
#include "tinyxml2.h"
#include "JetScapeXML.h"
#include "JetScapeWriter.h"

namespace Jetscape {


class SoftParticlization: public JetScapeModuleBase {
 private:

 public:
    SoftParticlization();
    ~SoftParticlization();

    virtual void Init();
    virtual void Exec();
    virtual void Clear();
    
    std::vector<std::vector<shared_ptr<Hadron>>> Hadron_list_;

    tinyxml2::XMLElement *xml_;
};

} // end namespace Jetscape


#endif  // SOFTPARTICLIZATION_H_
