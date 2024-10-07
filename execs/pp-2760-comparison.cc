/*
Code to compile the outputs of the analysis-spectra code for different Q0
and virtuality factors.
*/

//C++ header
#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"
#include <filesystem>

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
    vector<string> directories = getComparisonDirs(argc, argv); int maxdir = directories.size();
    cout << "Got point directories" << endl;
    chdir("/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/build");

    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //making said dat filesmakeObsPred(directories, input, "charged-xp", "hadrons", true);
    makeObsPred(directories, input, "HadronSpectraPredictionSoft", "soft hadrons");
    makeObsPred(directories, input, "HadronSpectraPredictionHard", "hard hadrons");
    makeObsPred(directories, input, "PionSpectraPredictionSoft", "rough soft pions");
    makeObsPred(directories, input, "PionSpectraPredictionHard", "rough hard pions");
    makeObsPred(directories, input, "KaonSpectraPredictionSoft", "rough soft kaons");
    makeObsPred(directories, input, "KaonSpectraPredictionHard", "rough hard kaons");
    makeObsPred(directories, input, "ProtonSpectraPredictionSoft", "rough soft protons");
    makeObsPred(directories, input, "ProtonSpectraPredictionHard", "rough hard protons");
    makeObsPred(directories, input, "JetSpectraPredictionR3", "jet radius 0.3");
    makeObsPred(directories, input, "JetSpectraPredictionR2", "jet radius 0.2");
    makeObsPred(directories, input, "JetSpectraPredictionR4", "jet radius 0.4");

    //exit behavior
    return 0;
}
