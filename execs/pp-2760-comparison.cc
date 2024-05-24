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

    //declaring vectors to hold results
    vector<vector<double>> jetPredicts, jetPredicts2, jetPredicts3;
    vector<vector<double>> hadronPredicts, pionPredicts, kaonPredicts, protonPredicts;
 
    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < maxdir; i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistHadrons = (TH1D*)file.Get("hadrons"); cout << "Got hadrons. ";
        TH1D* tempHistPions = (TH1D*)file.Get("smooth pions"); cout << "Got pions. ";
        TH1D* tempHistKaons = (TH1D*)file.Get("smooth kaons"); cout << "Got kaons. ";
        TH1D* tempHistProtons = (TH1D*)file.Get("smooth protons"); cout << "Got protons. ";
        TH1D* tempHistJets = (TH1D*)file.Get("jet radius 0.3"); cout << "Got jets 1. ";
        TH1D* tempHistJets2 = (TH1D*)file.Get("jet radius 0.2"); cout << "Got jets 2. ";
        TH1D* tempHistJets3 = (TH1D*)file.Get("jet radius 0.4"); cout << "Got jets 3." << endl;

        //adding to jets graph
        tempHistJets->Scale(1e6);
        tempHistJets2->Scale(1e6);
        tempHistJets3->Scale(1e6);

        //dat file data for Bayes analysis
        vector<double> jettemp, jettemp2, jettemp3, hadrontemp, piontemp, kaontemp, protontemp;
        for(int j = 1; j <= tempHistJets->GetNbinsX(); j++) jettemp.push_back(tempHistJets->GetBinContent(j));
        for(int j = 1; j <= tempHistJets2->GetNbinsX(); j++) jettemp2.push_back(tempHistJets2->GetBinContent(j));
        for(int j = 1; j <= tempHistJets3->GetNbinsX(); j++) jettemp3.push_back(tempHistJets3->GetBinContent(j));
        for(int j = 1; j <= tempHistHadrons->GetNbinsX(); j++) hadrontemp.push_back(tempHistHadrons->GetBinContent(j));
        for(int j = 1; j <= tempHistPions->GetNbinsX(); j++) piontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistKaons->GetNbinsX(); j++) kaontemp.push_back(tempHistKaons->GetBinContent(j));
        for(int j = 1; j <= tempHistProtons->GetNbinsX(); j++) protontemp.push_back(tempHistProtons->GetBinContent(j));
        jetPredicts.push_back(jettemp);
        jetPredicts2.push_back(jettemp2);
        jetPredicts3.push_back(jettemp3);
        hadronPredicts.push_back(hadrontemp);
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        
        file.Close();
    }

    //making said dat files
    chdir(input.c_str());
    char tmp[256]; getcwd(tmp, 256); cout << "Current working directory: " << tmp << endl;
    makeDatFile(jetPredicts, "JetSpectraPredictionR3", "# Version 1.0\n# Jet Spectra radius 0.3 for parameters.txt");
    makeDatFile(jetPredicts2, "JetSpectraPredictionR2", "# Version 1.0\n# Jet Spectra radius 0.2 for parameters.txt");
    makeDatFile(jetPredicts3, "JetSpectraPredictionR4", "# Version 1.0\n# Jet Spectra radius 0.4 for parameters.txt");
    makeDatFile(hadronPredicts, "HadronSpectraPrediction", "# Version 1.0\n# Hadron Spectra for parameters.txt");
    makeDatFile(pionPredicts, "PionSpectraPrediction", "# Version 1.0\n# Pion Spectra for parameters.txt");
    makeDatFile(kaonPredicts, "KaonSpectraPrediction", "# Version 1.0\n# Kaon Spectra for parameters.txt");
    makeDatFile(protonPredicts, "ProtonSpectraPrediction", "# Version 1.0\n# Proton Spectra for parameters.txt");

    //exit behavior
    return 0;
}
