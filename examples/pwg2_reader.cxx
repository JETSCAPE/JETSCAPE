/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
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

#include <TH1.h>
#include <TFile.h>

using namespace std;
using namespace Jetscape;

using namespace fjcore;

// -------------------------------------

// Forward declaration
void AnalyzeGraph(shared_ptr<PartonShower> mS);
ostream & operator<<(ostream & ostr, const PseudoJet & jet);

// -------------------------------------

// Create a pdf of the shower graph:
// Use with graphviz (on Mac: brew install graphviz --with-app)
// in shell: dot GVfile.gv -Tpdf -o outputPDF.pdf
// [or you can also use the GraphViz app for Mac Os X (in the "cellar" of homebrew)]

// -------------------------------------

int main(int argc, char** argv)
{

  // Parse arguments
  string inname  = "test_out.dat.gz";
  string outname = "test_out.root";
  string filetype = ".dat.gz";
  std::size_t pos=0;
  switch ( argc ){
    break;
  case 3:
    inname = argv[1];
    outname = argv[2];
    break;
  case 2:
    inname = argv[1];
    outname = inname;
    pos = outname.find( filetype);
    if ( pos==std::string::npos ){
      throw std::runtime_error ("Please use a filename ending in .dat.gz for an automatically generated output name");
      return -1;
    }
    outname.replace(pos, filetype.length(), ".root");
    break;
  case 1:
  case 0:
    break;    
  default:
    cout << "Usage: readerTestWithRoot [inputname] [outputname]" << endl;
    return -1;
  } 
    
  JetScapeLogger::Instance()->SetInfo(false);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  JetScapeLogger::Instance()->SetVerboseLevel(0);
  
  cout<<endl;

  // VERY crude way to pull event information that
  // resides at the foot of the file in a comment section
  // (for technical reasons);
  igzstream inFile;
  inFile.open(inname.c_str());
  string line;
  int nAccepted=-1;
  double sigmaGen=-1;
  string cs="# ";
  string tofind;
  while (getline(inFile,line)) {
    if ( line.compare(0, cs.length(), cs) != 0) continue;

    tofind = "# nAccepted = ";
    pos = line.find( tofind );
    if ( pos!=std::string::npos )      nAccepted=atoi( string ( line.begin()+pos+tofind.length(), line.end()).c_str());

    tofind = "# sigmaGen  = ";
    pos = line.find( tofind );
    if ( pos!=std::string::npos )      sigmaGen=atof( string ( line.begin()+pos+tofind.length(), line.end()).c_str());
  }
  if ( sigmaGen<0 || nAccepted<0 ){
    throw std::runtime_error("Couldn't find sigmaGen and nAccepted in the footer.");
    return -1;
  }

  INFO << "Using sigmaGen = " << sigmaGen << " and nAccepted = " << nAccepted;

  double weight = sigmaGen / nAccepted;
  // weight =1;
  
  // Simple jetfinding ...
  float R=0.4;
  JetDefinition jet_def(antikt_algorithm, R);

  // Constituent selectors
  // ---------------------
  float max_track_rap = 2;
  Selector select_track_rap = SelectorAbsRapMax(max_track_rap);
  Selector select_pt        = SelectorPtMin( 0.2 );
  Selector select_cons      = select_track_rap * select_pt; 

  // jet selectors
  Selector select_jet_eta = SelectorAbsEtaMax( max_track_rap-R );
  Selector select_jet_pt  = SelectorPtRange( 10,1000 );
  Selector select_jet     = select_jet_eta * select_jet_pt;
  

  // dijet cuts
  float ptlead = 40;
  float ptsublead = 20;
  
  auto reader=make_shared<JetScapeReaderAsciiGZ>(inname);
  
  // Open output root file
  TFile* fout = new TFile( outname.data(), "RECREATE");

  // turn on error bars by default
  TH1::SetDefaultSumw2();
  
  // Histograms
  TH1D* partonpt = new TH1D ("partonpt","Final Parton p_{T};p_{T} [GeV/c]",120, 0, 120);
  TH1D* partoneta = new TH1D ("partoneta","Final Parton #eta; #eta",100, -5, 5);
  TH1D* hardpartonpt = new TH1D ("hardpartonpt","Hard Parton p_{T};p_{T} [GeV/c]",80, 0, 40);
  TH1D* hardpartoneta = new TH1D ("hardpartoneta","Hard Parton #eta; #eta",100, -5, 5);

  TH1D* leadjetpt = new TH1D ("leadjetpt","Leading Jet p_{T};p_{T}^{jet} [GeV/c]",120, 0, 120);
  TH1D* subleadjetpt = new TH1D ("subleadjetpt","Sub-Leading Jet p_{T};p_{T}^{jet} [GeV/c]",120, 0, 120);

  TH1D* aj = new TH1D ("aj","A_{J};", 40, -1, 1 );

  TH1D* jetpt = new TH1D ("jetpt","Jet p_{T};p_{T}^{jet} [GeV/c]",120,0,120);
  TH1D* jeteta = new TH1D ("jeteta","Jet #eta;#eta^{jet} [GeV/c]",100, -3, 3);
    
  // Histogram for jet fragmentation function
  TH1D* jetz = new TH1D ("jetz","Jet z;p_{T}^{constituent} / p_{T}^{jet}",100, 0, 1);
    
  TH1D* dR = new TH1D ("dR","#DeltaR;counts",50, 0, 1);
  
  // reads in multiple events and multiple shower per event
  vector<shared_ptr<PartonShower>> mShowers;

  long event = 0; 
  while (!reader->Finished()){
    reader->Next();
    if ( !(event%1000) ) cout << "Working on " << event << endl;
    event++;
      
    // INFO<<"Analyze current event = "<<reader->GetCurrentEvent();
    mShowers=reader->GetPartonShowers();     

    // Collect jets from all showers
    vector<PseudoJet> jets;
    for (int i=0;i<mShowers.size();i++){	
      // Actual jetfinding
      vector<PseudoJet> particles = select_cons(mShowers[i]->GetFinalPartonsForFastJet());
      ClusterSequence cs(particles, jet_def);
	
      // mShowers[i]->PrintVertices();
      // mShowers[i]->PrintPartons();	
	
      // inital hard partons
      hardpartonpt->Fill ( mShowers[i]->GetPartonAt(0)->pt(), weight );
      hardpartoneta->Fill ( mShowers[i]->GetPartonAt(0)->eta(), weight );
	
      // final particles
      for ( auto &p : mShowers[i]->GetFinalPartonsForFastJet() ){
	partonpt->Fill ( p.pt(), weight );
	partoneta->Fill ( p.eta(), weight );
      }

      // jets
      for (auto j : select_jet( cs.inclusive_jets()) ){
	// additional cuts?
	jets.push_back(j);

	// save some constituent information
	for ( auto &consts : j.constituents() ){
	  // if ( consts.pt() > 40 ) {
	  //   cout << "---" << endl;
	  //   cout << consts << endl;
	  //   // cout << " --> " << *(mShowers[i]->GetPartonAt(0)) << endl;
	  //   cout << " --> shower initiator at pt = " << mShowers[i]->GetPartonAt(0)->pt()
	  // 	   << " eta = " << mShowers[i]->GetPartonAt(0)->eta()  << endl;
	  //   // for ( auto &parts : mShowers[i]->GetFinalPartons() ){
	  //   //   cout << " ----> " << *parts << endl;
	  //   //   cout << " ----> " << parts->pt() << endl;
	  //   // }
	      
	  // }
	  jetz->Fill( consts.pt() / j.pt(), weight );
	  dR->Fill( consts.delta_R(j), weight );
	}

      }


    }

    // analyze jets
    if ( jets.size()==0 ) continue;

    jets = sorted_by_pt( jets );
			   
    // misc information
    for (auto j : jets){
      // Fill histograms
      jetpt->Fill( j.pt(), weight );
      jeteta->Fill( j.eta(), weight );
	
    }

    // Dijets      
    if ( jets.size()<2 ) continue;
    PseudoJet lead = jets.at(0);
    PseudoJet sublead = jets.at(1);

    if ( lead.pt()<ptlead || sublead.pt()<ptsublead ) continue;

    // back to back?
    float dphi = lead.delta_phi_to( sublead ); // always in -pi,pi
    if ( dphi < -Jetscape::pi || dphi > Jetscape::pi ) {
      cerr << " dphi = " << dphi << endl;
      return -1;
    }

    if ( dphi > -Jetscape::pi+R && dphi < Jetscape::pi-R ) continue;

    // Fill histos
    leadjetpt->Fill( lead.pt(), weight );
    subleadjetpt->Fill( sublead.pt(), weight );
    aj->Fill( (lead.pt()-sublead.pt()) / (lead.pt()+sublead.pt()), weight);     
      
      
  }
  reader->Close();
  fout->Write();

  INFO_NICE<<"Finished!";
  INFO << "Read from " <<  inname;
  INFO << "Wrote to " <<  outname;
  cout << endl;
  return 0;

}

// -------------------------------------

void AnalyzeGraph(shared_ptr<PartonShower> mS)
{
  INFO<<"Some GTL graph/shower analysis/dfs search output:";
  
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

//----------------------------------------------------------------------
/// overloaded jet info output
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
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
