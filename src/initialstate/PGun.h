/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
//Parton Gun

#ifndef PGUN_H
#define PGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

class PGun: public HardProcess {
   
 private:
<<<<<<< HEAD
    double fixed_E0;
=======
    double fixed_pT;
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
    int parId;

 public:
    
    PGun();
     ~PGun();
     
     void InitTask();
     void Exec();
     
};

#endif  // PGUN_H
