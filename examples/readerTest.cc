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
  // commentend out so that you get the dot graph file for the first shower in the first event
  // (add in and the file gets overriden)
  //while (!reader->Finished())
    {
      reader->Next();

      cout<<"Analyze current event = "<<reader->GetCurrentEvent()<<endl;
      mShowers=reader->GetPartonShowers();     

      int finals = 0;
      for (int i=0;i<mShowers.size();i++)
	{
	  cout<<" Analyze parton shower = "<<i<<endl;
	 
	  mShowers[i]->PrintVertices();
	  mShowers[i]->PrintPartons();	

	  finals += mShowers[i]->GetFinalPartonsForFastJet().size();
	   
	  fjcore::ClusterSequence cs(mShowers[i]->GetFinalPartonsForFastJet(), jet_def);

	  vector<fjcore::PseudoJet> jets = fjcore::sorted_by_pt(cs.inclusive_jets(2));
	  cout<<endl;
	  cout<<jet_def.description()<<endl;
	  // Output of found jets ...
	  //cout<<endl;	 
	  for (int k=0;k<jets.size();k++)	    
	    cout<<"Anti-kT jet "<<k<<" : "<<jets[k]<<endl;
	  cout<<endl;
	  cout<<"Shower initiating parton : "<<*(mShowers[i]->GetPartonAt(0))<<endl;
	  cout<<endl;
	  
	  AnalyzeGraph(mShowers[i]);

	  if (i==0)
	    {
	      mShowers[i]->SaveAsGV("my_test.gv");
	      mShowers[i]->SaveAsGML("my_test.gml");
	      mShowers[i]->SaveAsGraphML("my_test.graphml");
	    }

	  // wait for 5s
	  //std::this_thread::sleep_for(std::chrono::milliseconds(5000));  
	}
      cout << " Found " << finals << " final state partons." << endl;
      auto hadrons = reader->GetHadrons();
      cout<<"Number of hadrons is: " << hadrons.size() << endl;
      
      fjcore::ClusterSequence hcs(reader->GetHadronsForFastJet(), jet_def);
      vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets(2));
      cout<<"AT HADRONIC LEVEL " << endl;
      for (int k=0;k<hjets.size();k++)	    
	cout<<"Anti-kT jet "<<k<<" : "<<hjets[k]<<endl;

      // for(unsigned int i=0; i<hadrons.size(); i++) {
      // 	cout<<"For Hadron Number "<<i<<" "<< hadrons[i].get()->e() << " "<< hadrons[i].get()->px()<< " "<< hadrons[i].get()->py() << " "<< hadrons[i].get()->pz()<< " "<< hadrons[i].get()->pt()<<  endl;
      // }
    }
    
    reader->Close(); 
}

// -------------------------------------

void AnalyzeGraph(shared_ptr<PartonShower> mS)
{
  INFO<<"Some GTL graph/shower analysis/dfs search output:";

  // quick and dirty ...
  dfs search;
  search.calc_comp_num(true);
  search.scan_whole_graph(true);
  search.start_node();// defaulted to first node ...
  search.run(*mS);

  cout<<endl;
  cout<<"DFS graph search feature from GTL:"<<endl;
  cout<<"Number of Nodes reached from node 0 = "<< search.number_of_reached_nodes () <<endl;
  cout<<"Node/Vertex ordering result from DFS:"<<endl;
  dfs::dfs_iterator itt2, endt2;
  for (itt2 = search.begin(), endt2=search.end(); itt2 !=endt2; ++itt2)
    {
      cout<<*itt2<<" ";//<<"="<<search.dfs_num(*itt2)<<" ";
    }
  cout<<endl;
  cout<<"Edge/Parton ordering result from DFS:"<<endl;
  dfs::tree_edges_iterator itt, endt;
  for (itt = search.tree_edges_begin(), endt=search.tree_edges_end(); itt !=endt; ++itt)
    {
        cout<<*itt;//<<endl;
    }
  cout<<endl;
  
  dfs::roots_iterator itt3, endt3;
  cout<<"List of root nodes found in graph/shower : ";
  for (itt3 = search.roots_begin(), endt3=search.roots_end(); itt3 !=endt3; ++itt3)
    {
      cout<<**itt3;
    }
  cout<<endl;
  cout<<endl;
}


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
