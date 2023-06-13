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
    TH1D* tempHistPions = (TH1D*)run.Get("raw pions"); cout << "Got pions. " << endl;

    //smoothing
    smoothBins(tempHistPions);
    //tempHistPions->Smooth();

    //writing out
    tempHistPions->Write("identified pions",TFile::kOverwrite);

    run.Close();
}