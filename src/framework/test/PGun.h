// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

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
