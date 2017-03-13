//Parton Gun Test ...

#ifndef PGUN_H
#define PGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"

class PGun: public HardProcess {
    // this is wrapper class for a simple brick
    // so that it can be used within the JETSCAPE framework
 private:
    double fixed_pT;

 public:
    
    PGun();
     ~PGun();
     
     void InitTask();
     void Exec();
     
};

#endif  // PGUN_H
