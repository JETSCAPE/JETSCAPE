/*
Smooths the identified hadrons graphs for an individual design point
./pp-smoothing /design/point/folder/
*/

//C++ header
#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"

#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"
#include "Pythia8/Pythia.h"

#include <GTL/dfs.h>
//ROOT headers
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TStyle.h"

#include "analysis.h"

using namespace std;
using namespace Jetscape;

int main(int argc, char* argv[]){
    // Create the ROOT application environment.
    chdir(argv[1]);
    TApplication theApp("hist", &argc, argv);

    //opening file and raw hists
    TFile run("root/totals.root","UPDATE");
    TH1D* tempHistPions = (TH1D*)run.Get("raw pions"); cout << "Got pions. ";
    TH1D* tempHistKaons = (TH1D*)run.Get("raw kaons"); cout << "Got kaons. ";
    TH1D* tempHistProtons = (TH1D*)run.Get("raw protons"); cout << "Got protons. " << endl;

    //smoothing
    smoothBins(tempHistPions);
    smoothBins(tempHistKaons);
    smoothBins(tempHistProtons);
    tempHistPions->Smooth();
    tempHistKaons->Smooth();
    tempHistProtons->Smooth();

    //writing out
    tempHistPions->Write("identified pions",TFile::kOverwrite);
    tempHistKaons->Write("identified kaons",TFile::kOverwrite);
    tempHistProtons->Write("identified protons",TFile::kOverwrite);

    //fixing jet cross sections
    /*TH1D* jet2 = (TH1D*)run.Get("jet radius 0.2"); cout << "Got R = 2. ";
    TH1D* jet3 = (TH1D*)run.Get("jet radius 0.3"); cout << "Got R = 3. ";
    TH1D* jet4 = (TH1D*)run.Get("jet radius 0.4"); cout << "Got R = 4. " << endl;

    jet2->Scale(62.8/43.1);
    jet3->Scale(62.8/43.1);
    jet4->Scale(62.8/43.1);

    jet2->Write("jet radius 0.2",TFile::kOverwrite);
    jet3->Write("jet radius 0.3",TFile::kOverwrite);
    jet4->Write("jet radius 0.4",TFile::kOverwrite);*/

    run.Close();
}