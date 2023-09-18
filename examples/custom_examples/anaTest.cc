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

// =======================================
// Put full FastJet before to avaoid binary/include errors!!!!
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"

#include "fastjet/Selector.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

// =======================================

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"

// =======================================

#include <Riostream.h>
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

using namespace std;
using namespace Jetscape;

const double etaMax  = 1.;

// -------------------------------------

// Forward declaration
void Show();
ostream & operator<<(ostream & ostr, const fastjet::PseudoJet & jet);

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

  string fName = "test_out";

  if (argc > 1)
    fName = argv[1];

  string outFile = (fName + ".root");
  //cout<<outFile<<endl;
  string inFile = (fName + ".dat.gz");
  //----------------------------------------------------------------------

  double z_cut = 0.1;
  double beta  = 0.0;

  fastjet::contrib::SoftDrop sd(beta, z_cut);
  cout << endl;
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  cout << "beta  = " << beta << endl;
  cout << "z_cut = " << z_cut << endl;

  //----------------------------------------------------------------------

  // first get some anti-kt jets
  double R = .4;
  double ptmin = 5.;

  // jet selectors
  fastjet::Selector sel_rap_jet = fastjet::SelectorAbsRapMax(etaMax-R);
  fastjet::Selector sel_match = fastjet::SelectorCircle(R);

  //----------------------------------------------------------------------

  TFile *f=new TFile(outFile.c_str(),"RECREATE");
  TH1D *hPt=new TH1D("hPt","",60,0,60);hPt->Sumw2();
  TH1D *hPtP=new TH1D("hPtP","",60,0,60);hPtP->Sumw2();
  TH1D *hM=new TH1D("hM","",80,0,20);hM->Sumw2();
  TH1D *hRg=new TH1D("hRg","",20,0,1);hRg->Sumw2();
  TH1D *hZg=new TH1D("hZg","",12,0,0.6);hZg->Sumw2();

  TH1D *hPtH=new TH1D("hPtH","",60,0,60);hPtH->Sumw2();
  TH1D *hMH=new TH1D("hMH","",80,0,20);hMH->Sumw2();
  TH1D *hRgH=new TH1D("hRgH","",20,0,1);hRgH->Sumw2();
  TH1D *hZgH=new TH1D("hZgH","",12,0,0.6);hZgH->Sumw2();

  //Do some dummy jetfinding ...
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  cout<<endl;
  cout<<jet_def.description()<<endl;
  cout<<endl;

  vector<shared_ptr<PartonShower>> mShowers;
  vector<fastjet::PseudoJet> vfinals;

  //auto reader=make_shared<JetScapeReaderAscii>("test_out_music.dat");
  auto reader=make_shared<JetScapeReaderAsciiGZ>(inFile);

  // reads in multiple events and multiple shower per event
  // commentend out so that you get the dot graph file for the first shower in the first event
  // (add in and the file gets overriden)

  int nJetsP = 0;
  int nJetsH = 0;
  int nEvents = 0;

  while (!reader->Finished())
    {
      reader->Next();

      //cout<<"Analyze current event = "<<reader->GetCurrentEvent()<<endl;
      mShowers=reader->GetPartonShowers();

      int finals = 0;
      for (int i=0;i<mShowers.size();i++)
      {
        //cout<<" Analyze parton shower = "<<i<<endl;
        //cout<<"Shower initiating parton : "<<*(mShowers[i]->GetPartonAt(0))<<endl;
        //cout<<endl;
        hPtP->Fill(mShowers[i]->GetPartonAt(0)->pt());

        //mShowers[i]->PrintVertices();
        //mShowers[i]->PrintPartons();

        finals += mShowers[i]->GetFinalPartonsForFastJet().size();
        auto mSfinal = mShowers[i]->GetFinalPartonsForFastJet();
        vfinals.insert(vfinals.end(),mSfinal.begin(), mSfinal.end());
      }

      //cout << " Found " << finals << " final state partons." << endl;
      //DEBUG:
      //cout<<vfinals.size()<<endl;

      fastjet::ClusterSequence cs(vfinals, jet_def);
      vector<fastjet::PseudoJet> jets = sel_rap_jet(fastjet::sorted_by_pt(cs.inclusive_jets(ptmin)));

      // Output of found jets ...
      //cout<<endl;
      /*
      sel_match.set_reference(jetDet);
      vector<fastjet::PseudoJet> jetPar_matched;
      jetPar_matched = fastjet::sorted_by_pt(sel_match(jets));
      */
      //cout<<"AT PARTONIC LEVEL " << endl;
      for (int k=0;k<jets.size();k++) {
        //cout<<"Anti-kT jet "<<k<<" : "<<jets[k]<<endl;

        double rg=sd(jets[k]).structure_of<fastjet::contrib::SoftDrop>().delta_R();
	      double zg=sd(jets[k]).structure_of<fastjet::contrib::SoftDrop>().symmetry();
        //cout<<rg<<" "<<zg<<endl;

        hPt->Fill(jets[k].pt());
        hM->Fill(jets[k].m());
        hRg->Fill(rg);
        hZg->Fill(zg);

        nJetsP++;

        //if (k>1) break;
      }

      //On Hadronic level ...
      auto hadrons = reader->GetHadrons();
      //cout<<"Number of hadrons is: " << hadrons.size() << endl;

      fastjet::ClusterSequence hcs(reader->GetHadronsForFastJet(), jet_def);
      vector<fastjet::PseudoJet> hjets = sel_rap_jet(fastjet::sorted_by_pt(hcs.inclusive_jets(ptmin)));

      //cout<<"AT HADRONIC LEVEL " << endl;
      for (int k=0;k<hjets.size();k++) {
        //cout<<"Anti-kT jet "<<k<<" : "<<hjets[k]<<endl;

        double rg=sd(hjets[k]).structure_of<fastjet::contrib::SoftDrop>().delta_R();
	      double zg=sd(hjets[k]).structure_of<fastjet::contrib::SoftDrop>().symmetry();
        //cout<<rg<<" "<<zg<<endl;

        hPtH->Fill(hjets[k].pt());
        hMH->Fill(hjets[k].m());
        hRgH->Fill(rg);
        hZgH->Fill(zg);

        nJetsH++;
      }

      //REMARK JP: Fill the histograms can be done smarter ...
      nEvents++;
      //cout<<endl;

      vfinals.clear();
      mShowers.clear();
    }

    f->Write();
    f->Close();

    reader->Close();

    cout<<"Analyzed #Events = "<<nEvents<<endl;
    cout<<"pTmin = "<<ptmin<<endl;
    cout<<"# Partonic Jets = "<<nJetsP<<endl;
    cout<<"# Hadronic Jets = "<<nJetsH<<endl;
    cout<<endl;
}

// -------------------------------------

void Show()
{
  ShowJetscapeBanner();
  INFO_NICE;
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"|  Ana Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}

//----------------------------------------------------------------------
/// overloaded jet info output

ostream & operator<<(ostream & ostr, const fastjet::PseudoJet & jet) {
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
