/*
Code to compile the outputs of the analysis-spectra code for different Q0
and virtuality factors.
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
    //list of directories to go over
    string input = argv[1]; 
    string pointsdir = input+"points/";
    vector<string> directories = getComparisonDirs(argc, argv);

    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //outputting results
    makeObsPred(directories, input, "charged-xp", "hadrons", true);
    makeObsPred(directories, input, "pion-xp", "pions", true);
    makeObsPred(directories, input, "kaon-xp", "kaons", true);
    makeObsPred(directories, input, "proton-xp", "protons", true);
    makeObsPred(directories, input, "jet", "jets", true);
    makeObsPred(directories, input, "dijet", "dijets", true);
    makeObsPred(directories, input, "mult", "multiplicity", true);
    return 0;
}
