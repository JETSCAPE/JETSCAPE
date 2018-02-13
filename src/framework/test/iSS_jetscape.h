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

 public:
    iSS_CF();
    ~iSS_CF();

    void Init();
    void Exec();
    void Clear();
};

#endif  // TEST_ISS_JETSCAPE_H_
