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

#include <GTL/dfs.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <thread>

#include "JetScapeBanner.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "PartonShower.h"
#include "fjcore.hh"
#include "gzstream.h"

using namespace std;
using namespace fjcore;

using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);

// -------------------------------------

int main(int argc, char** argv) {
  // JetScapeLogger::Instance()->SetInfo(false);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  // //SetVerboseLevel (9 a lot of additional debug output ...)
  // //If you want to suppress it: use SetVerboseLevle(0) or max
  // SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  // Whether to write the new header (ie. v2), including xsec info.
  // To enable, pass anything as the third argument to enable this option.
  // Default: disabled.
  bool writeHeaderV2 = false;
  if (argc > 3) {
    writeHeaderV2 = static_cast<bool>(argv[3]);
    std::cout << "NOTE: Writing header v2, and final cross section and error "
                 "at EOF.\n";
  }

  // The seperator between particles depends on the header.
  std::string particleSeperator = " ";
  if (!writeHeaderV2) {
    particleSeperator = "\t";
  }

  auto reader = make_shared<JetScapeReaderAscii>(argv[1]);
  std::ofstream dist_output(
      argv[2]);  // Format is SN, PID, E, Px, Py, Pz, Eta, Phi
  int SN = 0, TotalPartons = 0;
  while (!reader->Finished()) {
    reader->Next();

    // cout<<"Analyze current event: "<<reader->GetCurrentEvent()<<endl;
    auto mShowers = reader->GetPartonShowers();

    TotalPartons = 0;
    for (int i = 0; i < mShowers.size(); i++) {
      TotalPartons = TotalPartons + mShowers[i]->GetFinalPartons().size();
    }

    if (TotalPartons > 0) {
      ++SN;
      if (writeHeaderV2) {
        // NOTE: Needs consistent "\t" between all entries to simplify parsing
        // later.
        dist_output << "#"
                    << "\t"
                    << "Event\t" << SN << "\t"
                    << "weight\t" << reader->GetEventWeight() << "\t"
                    << "EPangle\t" << reader->GetEventPlaneAngle() << "\t"
                    << "N_partons\t" << TotalPartons << "\t"
                    << "|"  // As a delimiter
                    << "\t"
                    << "N"
                    << "\t"
                    << "pid"
                    << "\t"
                    << "status"
                    << "\t"
                    << "E"
                    << "\t"
                    << "Px"
                    << "\t"
                    << "Py"
                    << "\t"
                    << "Pz"
                    << "\n";
      } else {
        dist_output << "#"
                    << "\t" << reader->GetEventPlaneAngle() << "\t"
                    << "Event" << SN << "ID\t" << TotalPartons << "\t"
                    << "pstat-EPx"
                    << "\t"
                    << "Py"
                    << "\t"
                    << "Pz"
                    << "\t"
                    << "Eta"
                    << "\t"
                    << "Phi"
                    << "\n";
      }

      for (int i = 0; i < mShowers.size(); i++) {
        // cout<<" Analyze parton shower: "<<i<<endl;
        //  Let's create a file
        for (int ipart = 0; ipart < mShowers[i]->GetFinalPartons().size();
             ++ipart) {
          Parton p = *mShowers[i]->GetFinalPartons().at(ipart);
          //            if(abs(p.pid())!=5) continue;

          dist_output << ipart << particleSeperator << p.pid()
                      << particleSeperator << p.pstat() << particleSeperator
                      << p.e() << particleSeperator << p.px()
                      << particleSeperator << p.py() << particleSeperator
                      << p.pz();

          // v2 drops eta and phi, so only include it for v1
          if (!writeHeaderV2) {
            dist_output << particleSeperator << p.eta() << particleSeperator
                        << p.phi();
          }

          // Finish up
          dist_output << "\n";
        }
      }
    }
  }
  // Write the final cross section and error if requested by using header v2
  if (writeHeaderV2) {
    // NOTE: Needs consistent "\t" between all entries to simplify parsing
    // later.
    dist_output << "#"
                << "\t"
                << "sigmaGen\t" << reader->GetSigmaGen() << "\t"
                << "sigmaErr\t" << reader->GetSigmaErr() << "\n";
  }
  reader->Close();
}
