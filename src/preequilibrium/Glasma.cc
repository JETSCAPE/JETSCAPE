/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "Glasma.h"

#include <fstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>

#include "JetScapeConstants.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

using Jetscape::hbarC;

// Register the module with the base class
RegisterJetScapeModule<Glasma> Glasma::reg("Glasma");

Glasma::Glasma() {
  preequilibrium_status_ = NOT_STARTED;
  SetId("Glasma");
}

void Glasma::InitializePreequilibrium() {
  preequilibrium_status_ = INIT;
  preequilibrium_tau_0_ =
      GetXMLElementDouble({"Preequilibrium", "Glasma", "tau0"});
  preequilibrium_tau_max_ =
      GetXMLElementDouble({"Preequilibrium", "Glasma", "taus"});
  dtau_ = GetXMLElementDouble({"Preequilibrium", "Glasma", "dtau"});
}

void Glasma::EvolvePreequilibrium() {
  VERBOSE(2) << "Initialize density profiles in Glasma ...";
  std::string IPGlasmaFileName =
      GetXMLElementText({"Preequilibrium", "Glasma", "input_filename_glasma"});
  VERBOSE(2) << "Read in IP-Glasma Tmunu ...";
  std::ifstream IPGFile(IPGlasmaFileName.c_str());
  if (!IPGFile.good()) {
    Jetscape::JSWARN << "Can not open " << IPGlasmaFileName;
    exit(1);
  }
  std::string tempString;
  std::getline(IPGFile, tempString);
  double dummy;
  IPGFile >> dummy;
  while (!IPGFile.eof()) {
    double e_local;
    double u[4];
    double pi[10];
    IPGFile >> dummy >> dummy;
    IPGFile >> e_local >> u[0] >> u[1] >> u[2] >> u[3];
    e_.push_back(e_local);
    P_.push_back(e_local / 3.);
    utau_.push_back(u[0]);
    ux_.push_back(u[1]);
    uy_.push_back(u[2]);
    ueta_.push_back(u[3]);
    for (int i = 0; i < 10; i++) {
      IPGFile >> pi[i];
      pi[i] *= hbarC;
    }
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
    IPGFile >> dummy;
  }
  preequilibrium_status_ = DONE;
  IPGFile.close();
}
