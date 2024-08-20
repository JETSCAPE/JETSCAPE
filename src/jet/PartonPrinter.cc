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

#include "PartonPrinter.h"

#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include <iostream>
#include <string>

#include "GTL/graph.h"
#include "JetClass.h"
#include "JetScapeLogger.h"

using namespace std;

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<PartonPrinter> PartonPrinter::reg("PartonPrinter");

PartonPrinter::PartonPrinter() {}

PartonPrinter::~PartonPrinter() {}

void PartonPrinter::Init() {
  this->SetId("PartonPrinter");
  JSINFO << "Initialize PartonPrinter ...";
  string filename = GetXMLElementText({"PartonPrinter", "FileName"});

  dist_output.open(filename, std::ios::app);
}

void PartonPrinter::Exec() {
  VERBOSE(2) << "Run PartonPrinter: print shower from event # "
             << GetCurrentEvent() << " ...";
}

void PartonPrinter::GetFinalPartons(
    shared_ptr<PartonShower>
        pShower /*, vector<shared_ptr<Parton>>& fPartons*/) {
  if (pShower) {
    // vector<shared_ptr<Parton>> vPin;
    for (unsigned int ipart = 0;
         ipart < pShower.get()->GetFinalPartons().size(); ++ipart) {
      // fPartons.push_back( pShower.get()->GetFinalPartons().at(ipart));
      if (std::abs(pShower.get()->GetFinalPartons().at(ipart)->pid()) == 21) {
        dist_output << ipart << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->pid() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->e() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->px() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->py() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->pz() << endl;
      }

      // vPin.push_back( pShower.get()->GetFinalPartons().at(ipart));
    }
    // this->pFinals.push_back(vPin);
  }
}

void PartonPrinter::GetFinalPartons2(shared_ptr<PartonShower> pShower) {
  if (pShower) {
    /*
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size();
    ++ipart)
    {
      this->pFinals.push_back( pShower.get()->GetFinalPartons());
      //cout << "############### FINAL PARTON IN THE VECTOR NUMBER : " << ipart
    << " is " << this->pFinals.at(ipart) << "###################\n";
    }
*/
    // this->pFinals.clear();
    this->pFinals.push_back(pShower.get()->GetFinalPartons());
  }
}

void PartonPrinter::Clear() {
  // dist_output << " *********************************** " << endl;
  // dist_output.close();
  this->pFinals.clear();
}

void PartonPrinter::GetPartonsAtTime(shared_ptr<PartonShower> pShower,
                                     vector<shared_ptr<Parton>> &fPartons,
                                     double time) {}

}  // end namespace Jetscape
