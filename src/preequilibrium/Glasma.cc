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

#include <stdio.h>
#include <sys/stat.h>

#include <string>
#include <fstream>

#include "JetScapeLogger.h"
#include "Glasma.h"

// Register the module with the base class
RegisterJetScapeModule<Glasma> Glasma::reg("Glasma");

Glasma::Glasma() {
    preequilibrium_status_ = NOT_STARTED;
    SetId("Glasma");
}

void Glasma::EvolvePreequilibrium() {
    VERBOSE(2) << "Initialize density profiles in Glasma ...";
    std::string IPGlasmaFileName = "epsilon-u-hydro-t0.6-0.dat";
    VERBOSE(2) << "Read in IPGlasma Tmunu ...";
    std::ifstream IPGFile(IPGlasmaFileName.c_str());
    std::string tempString;
    std::getline(IPGFile, tempString);
    double dummy;
    IPGFile >> dummy >> dummy >> dummy;
    while (!IPGFile.eof()) {
        double e_local;
        double u[4];
        double pi[10];
        IPGFile >> e_local >> u[0] >> u[1] >> u[2] >> u[3];
        e_.push_back(e_local);
        P_.push_back(e_local/3.);
        utau_.push_back(u[0]);
        ux_.push_back(u[1]);
        uy_.push_back(u[2]);
        ueta_.push_back(u[3]);
        for (int i = 0; i < 10; i++) IPGFile >> pi[i];
        pi00_.push_back(pi[0]);
        pi01_.push_back(pi[1]);
        pi02_.push_back(pi[2]);
        pi03_.push_back(pi[3]);
        pi11_.push_back(pi[4]);
        pi12_.push_back(pi[5]);
        pi13_.push_back(pi[6]);
        pi22_.push_back(pi[7]);
        pi23_.push_back(pi[8]);
        pi33_.push_back(pi[9]);
        bulk_Pi_.push_back(0.);
        IPGFile >> dummy >> dummy >> dummy;
    }
    preequilibrium_status_ = DONE;
}
