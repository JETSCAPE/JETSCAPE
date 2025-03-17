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
    vector<string> directories = getComparisonDirs(argc, argv); int maxdir = directories.size();
    cout << "Got point directories" << endl;
    chdir("/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/build");

    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //making said dat files
    makeObsPred(directories, input, "PionSpectraPredictionSoft", "identified pions");
    makeObsPred(directories, input, "PionSpectraPredictionHard", "identified hard pions");
    makeObsPred(directories, input, "KaonSpectraPredictionSoft", "identified kaons");
    makeObsPred(directories, input, "ProtonSpectraPredictionSoft", "identified protons");
    makeObsPred(directories, input, "JetSpectraPrediction", "jets");

    //exit behavior
    return 0;
}
