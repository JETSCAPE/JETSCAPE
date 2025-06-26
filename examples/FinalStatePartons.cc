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

  // Whether to write a particular header version (eg. v2), including xsec info.
  // To enable, pass the desired version (just the number) as the third argument.
  // Default: v1
  unsigned int headerVersion = 1;
  if (argc > 3) {
    headerVersion = std::atoi(argv[3]);
    std::cout << "NOTE: Writing header v" << headerVersion << ", and final cross section and error at EOF.\n";
  }
  std::cout << "NOTE: Writing with output version v" << headerVersion << "\n";

  // The separator between particles depends on the header.
  std::string particleSeparator = " ";
  if (headerVersion == 1) {
    particleSeparator = "\t";
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
      if (headerVersion > 1) {
        // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
        dist_output << "#"
            << "\t" << "Event\t" << SN
            << "\t" << "weight\t" << reader->GetEventWeight()
            << "\t" << "EPangle\t" << reader->GetEventPlaneAngle()
            << "\t" << "N_partons\t" << TotalPartons;
        if (headerVersion == 3) {
          dist_output
              << "\t" << "vertex_x\t" << reader->GetVertexX()
              << "\t" << "vertex_y\t" << reader->GetVertexY()
              << "\t" << "vertex_z\t" << reader->GetVertexZ();
        }
        dist_output
            << "\t" << "|"  // As a delimiter
            << "\t" << "N"
            << "\t" << "pid"
            << "\t" << "status"
            << "\t" << "E"
            << "\t" << "Px"
            << "\t" << "Py"
            << "\t" << "Pz"
            << "\n";
      }
      else {
        dist_output << "#" << "\t"
            << reader->GetEventPlaneAngle() << "\t"
            << "Event"
            << SN << "ID\t"
            << TotalPartons << "\t"
            << "pstat-EPx"  << "\t"
            << "Py"  << "\t"
            << "Pz"  << "\t"
            << "Eta" <<  "\t"<< "Phi" << "\n";
      }

      for (int i = 0; i < mShowers.size(); i++) {
        // cout<<" Analyze parton shower: "<<i<<endl;
        //  Let's create a file
        for (int ipart = 0; ipart < mShowers[i]->GetFinalPartons().size();
             ++ipart) {
          Parton p = *mShowers[i]->GetFinalPartons().at(ipart);
          //            if(abs(p.pid())!=5) continue;

          dist_output << ipart << particleSeparator << p.pid()
                      << particleSeparator << p.pstat() << particleSeparator
                      << p.e() << particleSeparator << p.px()
                      << particleSeparator << p.py() << particleSeparator
                      << p.pz();

          // v2 drops eta and phi, so only include it for v1
          if (headerVersion == 1) {
            dist_output << particleSeparator << p.eta() << particleSeparator
                        << p.phi();
          }

          // Finish up
          dist_output << "\n";
        }
      }
    }
  }
  // Write the final cross section and error if requested by using header v2
  if (headerVersion > 1) {
    // NOTE: Needs consistent "\t" between all entries to simplify parsing
    // later.
    dist_output << "#"
                << "\t"
                << "sigmaGen\t" << reader->GetSigmaGen() << "\t"
                << "sigmaErr\t" << reader->GetSigmaErr() << "\n";
  }
  reader->Close();
}
