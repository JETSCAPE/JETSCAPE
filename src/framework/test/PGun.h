//Parton Gun Test ...

#ifndef PGUN_H
#define PGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"

class PGun: public HardProcess {
   
 private:
    double fixed_pT;

 public:
    
    PGun();
     ~PGun();
     
     void InitTask();
     void Exec();
     
};

#endif  // PGUN_H
