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
// Reader test (focus on graph)

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
//using namespace fjcore;

using namespace Jetscape;

// -------------------------------------

// Forward declaration
void Show();
void AnalyzeGraph(shared_ptr<PartonShower> mS);
ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);

// -------------------------------------

// Create a pdf of the shower graph:
// Use with graphviz (on Mac: brew install graphviz --with-app)
// in shell: dot GVfile.gv -Tpdf -o outputPDF.pdf
// [or you can also use the GraphViz app for Mac Os X (in the "cellar" of homebrew)]

// -------------------------------------

int main(int argc, char** argv)
{
  std::fstream fs;
  fs.open ("hqjet_pT.dat", std::fstream::out | std::fstream::app);
  
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);
  
  cout<<endl;
  Show();

  //Do some dummy jetfinding ...
  fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, 0.7);
  
  vector<shared_ptr<PartonShower>> mShowers;
  
  //Directly with template: provide the relevant stream
  //auto reader=make_shared<JetScapeReader<ifstream>>("test_out.dat");
  //auto reader=make_shared<JetScapeReader<igzstream>>("test_out.dat.gz");
  
  // Hide Template (see class declarations in reader/JetScapeReader.h) ...
  auto reader=make_shared<JetScapeReaderAscii>("test_out.dat");
  //auto reader=make_shared<JetScapeReaderAsciiGZ>("test_out.dat.gz");

  // reads in multiple events and multiple shower per event
  // commented out so that you get the dot graph file for the first shower in the first event
  // (add in and the file gets overriden)
  while (!reader->Finished())
    {
      reader->Next();

      cout<<"-----------------------------------------------------------"<<endl;
      cout<<"Analyze current event = "<<reader->GetCurrentEvent()<<endl;
      mShowers=reader->GetPartonShowers();     

      int finals = 0;
      vector<fjcore::PseudoJet> finalpartons;
      for (int i=0;i<mShowers.size();i++)
      {
        vector<fjcore::PseudoJet> temp = mShowers[i]->GetFinalPartonsForFastJet();
        finals += temp.size();
        for(int j=0; j<temp.size(); j++)
        {
	  finalpartons.push_back(temp[j]);  
        }
      }

      vector<string> heavyjets;

      fjcore::ClusterSequence cs(finalpartons, jet_def);
      vector<fjcore::PseudoJet> jets = fjcore::sorted_by_pt(cs.inclusive_jets(2));
      cout<<endl;
      cout<<jet_def.description()<<endl;
      cout<<endl;

      fs<<"event:"<<endl;
      for (int k=0;k<jets.size();k++)	    
      {
        cout<<"Anti-kT jet "<<k<<" : "<<jets[k]<<endl;
        vector<fjcore::PseudoJet> constituents = jets[k].constituents();
        for(int p=0; p<constituents.size(); p++)
        {
          int pid = constituents[p].user_info<Parton>().pid();
          cout<<"    pid: "<<pid<<" four vec: " << constituents[p] << " hq_channel: " << to_string(constituents[p].user_info<Parton>().hq_channel()) << endl;
	  
          if(abs(pid)==4||abs(pid)==5)
          {
            string heavyjet =to_string(pid)+ " " + to_string(jets[k].pt()) + " " + to_string(jets[k].m()) + " " + to_string(jets[k].rap()) + " " + to_string(jets[k].phi());
            heavyjet += " " + to_string(constituents[p].user_info<Parton>().hq_channel()) + " " + to_string(constituents[p].user_info<Parton>().hq_mother_id());
            heavyjets.push_back(heavyjet);
	    //print hqjet_info to a new file
	    fs<<heavyjet<<endl;
	  }
	  
        }
	cout<<endl;
      }
      
      cout << " Found " << finals << " final state partons." << endl;
      cout<<"-----------------------------------------------------------"<<endl;
      for(int i=0; i<heavyjets.size(); i++)
      {
         cout << heavyjets[i] << endl;
      }
      cout<<"-----------------------------------------------------------"<<endl;
    }
    
    reader->Close();
    fs.close();
}

// -------------------------------------



// -------------------------------------

void Show()
{
  ShowJetscapeBanner();
  INFO_NICE;
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"| Reader Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}

//----------------------------------------------------------------------
/// overloaded jet info output

ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi();
  }
  return ostr;
}


//----------------------------------------------------------------------
