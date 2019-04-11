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
    double fixed_pT;

    // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<PGun> reg;

 public:
    
    PGun();
     ~PGun();
     
     void InitTask();
     void Exec();
     
};

#endif  // PGUN_H
