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

    //data arrays
    vector<vector<double>> specPredicts;
    vector<vector<double>> thrustPredicts;
    vector<vector<double>> multiplicityPredicts;
    vector<vector<double>> pionPredicts;
    vector<vector<double>> kaonPredicts;
    vector<vector<double>> protonPredicts;
    vector<vector<double>> jetPredicts;
    vector<vector<double>> dijetPredicts;

    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < directories.size(); i++){
        //temp vectors
        vector<double> spectemp, thrusttemp, multtemp, piontemp, kaontemp, protontemp, jettemp, dijettemp;

        //reading graph
        string filename = directories[i] + "/totals.root";
        cout << directories[i] << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistHadronSpecP = (TH1D*)file.Get("hadrons;1"); cout << "Got hadrons. ";
        TH1D* tempHistThrustP = (TH1D*)file.Get("thrust;1"); cout << "Got thrust. ";
        TH1D* tempHistMultiplicityP = (TH1D*)file.Get("multiplicity;1"); cout << "Got multiplicity. ";
        TH1D* tempHistPionsP = (TH1D*)file.Get("pions;1"); cout << "Got pions. ";
        TH1D* tempHistKaonsP = (TH1D*)file.Get("kaons;1"); cout << "Got kaons. ";
        TH1D* tempHistProtonsP = (TH1D*)file.Get("protons;1"); cout << "Got protons. ";
        TH1D* tempHistJetsP = (TH1D*)file.Get("jets;1"); cout << "Got jets. ";
        TH1D* tempHistDiJetsP = (TH1D*)file.Get("dijets;1"); cout << "Got dijets. " << endl;

        //adding to temp vectors, with conditions for skipped bins in identified hadrons
        for(int j = 0; j < tempHistHadronSpecP->GetNbinsX(); j++) spectemp.push_back(tempHistHadronSpecP->GetBinContent(j+1));
        for(int j = 0; j < tempHistThrustP->GetNbinsX(); j++) thrusttemp.push_back(tempHistThrustP->GetBinContent(j+1));
        for(int j = 0; j < tempHistMultiplicityP->GetNbinsX(); j++) multtemp.push_back(tempHistMultiplicityP->GetBinContent(j+1));
        for(int j = 0; j < tempHistPionsP->GetNbinsX(); j++) if(tempHistPionsP->GetBinContent(j+1)>0) piontemp.push_back(tempHistPionsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistKaonsP->GetNbinsX(); j++) if(tempHistKaonsP->GetBinContent(j+1)>0) kaontemp.push_back(tempHistKaonsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistProtonsP->GetNbinsX(); j++) if(tempHistProtonsP->GetBinContent(j+1)>0) protontemp.push_back(tempHistProtonsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistJetsP->GetNbinsX(); j++) if(tempHistJetsP->GetBinContent(j+1)>0) jettemp.push_back(tempHistJetsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistDiJetsP->GetNbinsX(); j++) if(tempHistDiJetsP->GetBinContent(j+1)>0) dijettemp.push_back(tempHistDiJetsP->GetBinContent(j+1));

        //adding to full vectors
        specPredicts.push_back(spectemp);
        thrustPredicts.push_back(thrusttemp);
        multiplicityPredicts.push_back(multtemp);
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        jetPredicts.push_back(jettemp);
        dijetPredicts.push_back(dijettemp);

        file.Close();
    }

    //outputting results
    chdir(input.c_str());
    makeDatFile(specPredicts, "charged-xp", "# Version 1.0\n# Spectra for parameters.txt");
    makeDatFile(thrustPredicts, "thrust", "# Version 1.0\n# Thrust for parameters.txt");
    makeDatFile(multiplicityPredicts, "mult", "# Version 1.0\n# Multiplicity for parameters.txt");
    makeDatFile(pionPredicts, "pion-xp", "# Version 1.0\n# Pions for parameters.txt");
    makeDatFile(kaonPredicts, "kaon-xp", "# Version 1.0\n# Kaons for parameters.txt");
    makeDatFile(protonPredicts, "proton-xp", "# Version 1.0\n# Protons for parameters.txt");
    makeDatFile(jetPredicts, "jet", "# Version 1.0\n# Jets for parameters.txt");
    makeDatFile(dijetPredicts, "dijet", "# Version 1.0\n# Dijets for parameters.txt");

    return 0;
}
