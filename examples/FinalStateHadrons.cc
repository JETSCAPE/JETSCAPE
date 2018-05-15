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
  // //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);
  
  auto reader=make_shared<JetScapeReaderAscii>("test_out.dat");
  std::ofstream dist_output ("JetscapeFinalStateHadrons.txt"); //Format is SN, PID, E, Px, Py, Pz, Eta, Phi
  vector<shared_ptr<Hadron>> hadrons;
  while (!reader->Finished())
    {
      reader->Next();
      
      // cout<<"Analyze current event: "<<reader->GetCurrentEvent()<<endl;

      //dist_output<<"Event "<< reader->GetCurrentEvent()+1<<endl;
      hadrons = reader->GetHadrons();
      cout<<"Number of hadrons is: " << hadrons.size() << endl;
      for(unsigned int i=0; i<hadrons.size(); i++)
	{
	  dist_output<<i<<" "<<hadrons[i].get()->pid()<<" "<<hadrons[i].get()->pstat()<<" "<< hadrons[i].get()->e() << " "<< hadrons[i].get()->px()<< " "<< hadrons[i].get()->py() << " "<< hadrons[i].get()->pz()<<  endl;
        }
    }
  
  reader->Close();
  
}
