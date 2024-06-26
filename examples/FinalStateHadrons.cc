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

#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
using namespace fjcore;

using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);


// -------------------------------------

int main(int argc, char** argv)
{
  // JetScapeLogger::Instance()->SetInfo(false);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  // //SetVerboseLevel (9 a lot of additional debug output ...)
  // //If you want to suppress it: use SetVerboseLevel(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  // Whether to write the new header (ie. v2), including xsec info.
  // To enable, pass anything as the third argument to enable this option.
  // Default: version 1.
  int headerVersion = 1;
  if (argc > 3) {
    headerVersion = static_cast<int>(argv[3]);
  }
  std::cout << "NOTE: Writing with output version v" << headerVersion << "\n";

  // The separator between particles _does not_ depend on the header for final state hadrons
  std::string particleSeparator = " ";

  auto reader = make_shared<JetScapeReaderAscii>(argv[1]);
  std::ofstream dist_output (argv[2]); //Format is SN, PID, E, Px, Py, Pz, Eta, Phi
  vector<shared_ptr<Hadron>> hadrons;
  int SN=0;
  while (!reader->Finished())
  {
    reader->Next();

    // cout<<"Analyze current event: "<<reader->GetCurrentEvent()<<endl;

    //dist_output<<"Event "<< reader->GetCurrentEvent()+1<<endl;
    hadrons = reader->GetHadrons();
    cout<<"Number of hadrons is: " << hadrons.size() << endl;

    if (hadrons.size() > 0)
    {
      ++SN;
      if (headerVersion == 2) {
        // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
        dist_output << "#"
            << "\t" << "Event\t" << SN
            << "\t" << "weight\t" << reader->GetEventWeight()
            << "\t" << "EPangle\t" << reader->GetEventPlaneAngle()
            << "\t" << "N_hadrons\t" << hadrons.size()
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
            << hadrons.size() << "\t"
            << "pstat-EPx"   << "\t"
            << "Py"  << "\t"
            << "Pz"  << "\t"
            << "Eta" <<  "\t"<< "Phi" << "\n";
      }

      for (unsigned int i=0; i<hadrons.size(); i++)
      {
        dist_output << i
            << particleSeparator << hadrons[i].get()->pid()
            << particleSeparator << hadrons[i].get()->pstat()
            << particleSeparator << hadrons[i].get()->e()
            << particleSeparator << hadrons[i].get()->px()
            << particleSeparator << hadrons[i].get()->py()
            << particleSeparator << hadrons[i].get()->pz();

        // v2 drops eta and phi, so only include it for v1
        if (headerVersion == 1) {
            dist_output << particleSeparator << hadrons[i].get()->eta()
                << particleSeparator << hadrons[i].get()->phi();
        }

        // Finish up
        dist_output << "\n";
      }
    }
  }
  // Write the final cross section and error if requested by using header v2
  if (headerVersion > 1) {
    // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
    dist_output << "#"
        << "\t" << "sigmaGen\t" << reader->GetSigmaGen()
        << "\t" << "sigmaErr\t" << reader->GetSigmaErr()
        << "\n";
  }
  reader->Close();
}
