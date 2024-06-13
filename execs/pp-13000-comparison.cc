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
    vector<string> directories = getComparisonDirs(argc, argv);
    int maxdir = directories.size();
    string input = argv[1]; 
    
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //declaring vectors to hold results
    vector<vector<double>> pionPredicts, kaonPredicts, protonPredicts, jet1Predicts, jet2Predicts, jet3Predicts;
 
    //grabbing graphs from each jetscape file
    for(int i = 0; i < directories.size(); i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistPions = (TH1D*)file.Get("identified pions"); cout << "Got pions. ";
        TH1D* tempHistKaons = (TH1D*)file.Get("identified kaons"); cout << "Got kaons. ";
        TH1D* tempHistProtons = (TH1D*)file.Get("identified protons"); cout << "Got protons. ";
        TH1D* tempHistJets1 = (TH1D*)file.Get("low y jets"); cout << "Got low y jets. ";
        TH1D* tempHistJets2 = (TH1D*)file.Get("mid y jets"); cout << "Got mid y jets. ";
        TH1D* tempHistJets3 = (TH1D*)file.Get("high y jets"); cout << "Got high y jets. " << endl;

        //dat file data for Bayes analysis
        vector<double> piontemp, kaontemp, protontemp, jettemp1, jettemp2, jettemp3;
        for(int j = 1; j <= tempHistPions->GetNbinsX(); j++) piontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistKaons->GetNbinsX(); j++) kaontemp.push_back(tempHistKaons->GetBinContent(j));
        for(int j = 1; j <= tempHistProtons->GetNbinsX(); j++) protontemp.push_back(tempHistProtons->GetBinContent(j));
        for(int j = 1; j <= tempHistJets1->GetNbinsX(); j++) jettemp1.push_back(tempHistJets1->GetBinContent(j));
        for(int j = 1; j <= tempHistJets2->GetNbinsX(); j++) jettemp2.push_back(tempHistJets2->GetBinContent(j));
        for(int j = 1; j <= tempHistJets3->GetNbinsX(); j++) jettemp3.push_back(tempHistJets3->GetBinContent(j));
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        jet1Predicts.push_back(jettemp1);
        jet2Predicts.push_back(jettemp2);
        jet3Predicts.push_back(jettemp3);
        
        file.Close();
    }

    //making said dat files
    chdir(input.c_str());
    char tmp[256]; getcwd(tmp, 256); cout << "Current working directory: " << tmp << endl;
    makeDatFile(pionPredicts, "PionSpectraPrediction", "# Version 1.0\n# Pion Spectra for parameters.txt");
    makeDatFile(kaonPredicts, "KaonSpectraPrediction", "# Version 1.0\n# Kaon Spectra for parameters.txt");
    makeDatFile(protonPredicts, "ProtonSpectraPrediction", "# Version 1.0\n# Proton Spectra for parameters.txt");
    makeDatFile(jet1Predicts, "LowYJetPrediction", "# Version 1.0\n# Low Y Jet Spectra for parameters.txt");
    makeDatFile(jet2Predicts, "MidYJetPrediction", "# Version 1.0\n# Mid Y Jet Spectra for parameters.txt");
    makeDatFile(jet3Predicts, "HighYJetPrediction", "# Version 1.0\n# High Y Jet Spectra for parameters.txt");

    //exit behavior
    return 0;
}
