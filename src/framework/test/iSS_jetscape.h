// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// Copyright [2018] <Chun Shen>
// -----------------------------------------

#ifndef TEST_ISS_JETSCAPE_H_
#define TEST_ISS_JETSCAPE_H_

#include "SoftParticlization.h"
#include "iSS.h"

using namespace Jetscape;

class iSS_CF: public SoftParticlization {
 private:
    tinyxml2::XMLElement *iSS_xml_;

    iSS *iSpectraSampler_ptr_;

 public:
    iSS_CF();
    ~iSS_CF();

    void InitTask();
    void Exec();
    void Clear();
    void WriteTask(weak_ptr<JetScapeWriter> w);

    void pass_hadron_list_to_JETSCAPE();
};

#endif  // TEST_ISS_JETSCAPE_H_
