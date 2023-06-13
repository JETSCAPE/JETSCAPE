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
    vector<string> directories = {};
    string input = argv[1]; 

    for(int i = 1; i < argc; i++){
        chdir(argv[i]);
        vector<string> tempdirs = get_directories(".");

        //removing the results dir
        int rmindex = 0;
        for(int j = 0; j < tempdirs.size(); j++){
            tempdirs[j].erase(0,2);
            if(tempdirs[j].find("QVir_Analysis") != string::npos) rmindex = j;
        }
        tempdirs.erase(tempdirs.begin() + rmindex); 

        vector<string> sorteddirs = doubleSort(tempdirs); //sorting for consistency
        for(int j = 0; j < sorteddirs.size(); j++){
            string temp = argv[i] + sorteddirs[j];
            sorteddirs[j] = temp;
        }
        directories.insert(directories.end(), sorteddirs.begin(), sorteddirs.end()); //inserting into the end of total vector
    }
    chdir("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/build");

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

    //dijet comparison
    TMultiGraph* jets = new TMultiGraph("Jet Spectra","Jet Spectra");
    TGraphErrors* jetscapejets[directories.size()];
    TFile jet_file("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/data/eeDiJet.root");
    TDirectory* jetdir = (TDirectory*)jet_file.Get("LeadingDiJetEnergy");
    TGraphErrors* jetData = (TGraphErrors*)jetdir->Get("Graph1D_y1");
    jetData->SetMarkerStyle(kFullDotLarge);
    jets->Add(jetData,"ap");

    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < directories.size(); i++){
        //temp vectors
        vector<double> spectemp, thrusttemp, multtemp, piontemp, kaontemp, protontemp, jettemp;

        //reading graph
        string filename = directories[i] + "/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistHadronSpecP = (TH1D*)file.Get("hadrons;1"); cout << "Got hadrons. ";
        TH1D* tempHistThrustP = (TH1D*)file.Get("thrust;1"); cout << "Got thrust. ";
        TH1D* tempHistMultiplicityP = (TH1D*)file.Get("multiplicity;1"); cout << "Got multiplicity. ";
        TH1D* tempHistPionsP = (TH1D*)file.Get("pions;1"); cout << "Got pions. ";
        TH1D* tempHistKaonsP = (TH1D*)file.Get("kaons;1"); cout << "Got kaons. ";
        TH1D* tempHistProtonsP = (TH1D*)file.Get("protons;1"); cout << "Got protons. ";
        TH1D* tempHistJetsP = (TH1D*)file.Get("jets;1"); cout << "Got jets. " << endl;

        //adding to temp vectors, with conditions for skipped bins in identified hadrons
        for(int j = 0; j < tempHistHadronSpecP->GetNbinsX(); j++) spectemp.push_back(tempHistHadronSpecP->GetBinContent(j+1));
        for(int j = 0; j < tempHistThrustP->GetNbinsX(); j++) thrusttemp.push_back(tempHistThrustP->GetBinContent(j+1));
        for(int j = 0; j < tempHistMultiplicityP->GetNbinsX(); j++) multtemp.push_back(tempHistMultiplicityP->GetBinContent(j+1));
        for(int j = 0; j < tempHistPionsP->GetNbinsX(); j++) if(tempHistPionsP->GetBinContent(j+1)>0) piontemp.push_back(tempHistPionsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistKaonsP->GetNbinsX(); j++) if(tempHistKaonsP->GetBinContent(j+1)>0) kaontemp.push_back(tempHistKaonsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistProtonsP->GetNbinsX(); j++) if(tempHistProtonsP->GetBinContent(j+1)>0) protontemp.push_back(tempHistProtonsP->GetBinContent(j+1));
        for(int j = 0; j < tempHistJetsP->GetNbinsX(); j++) if(tempHistJetsP->GetBinContent(j+1)>0) jettemp.push_back(tempHistJetsP->GetBinContent(j+1));

        //adding to full vectors
        specPredicts.push_back(spectemp);
        thrustPredicts.push_back(thrusttemp);
        multiplicityPredicts.push_back(multtemp);
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        jetPredicts.push_back(jettemp);

        //dijet plotting
        TH1D* tempHistDiJets = (TH1D*)file.Get("dijets;1");
        jetscapejets[i] = (TGraphErrors*)histToGraph(tempHistDiJets).Clone();
        jets->Add(jetscapejets[i],"lX");

        file.Close();
    }

    //outputting results
    chdir(input.c_str());
    makeDatFile(specPredicts, "SpectrumPrediction", "# Version 1.0\n# Spectra for parameters.txt");
    makeDatFile(thrustPredicts, "thrustPrediction", "# Version 1.0\n# Thrust for parameters.txt");
    makeDatFile(multiplicityPredicts, "multiplicityPrediction", "# Version 1.0\n# Multiplicity for parameters.txt");
    makeDatFile(pionPredicts, "pionPrediction", "# Version 1.0\n# Pions for parameters.txt");
    makeDatFile(kaonPredicts, "kaonPrediction", "# Version 1.0\n# Kaons for parameters.txt");
    makeDatFile(protonPredicts, "protonPrediction", "# Version 1.0\n# Protons for parameters.txt");
    makeDatFile(jetPredicts, "eejetPrediction", "# Version 1.0\n# Jets for parameters.txt");
    
    //dijet drawing
    TCanvas* cJet = new TCanvas("c2","c2",1400,1200);
    cJet->cd();
    cJet->SetLeftMargin(0.15);
    cJet->SetLogy();
    jets->Draw();
    cJet->Print("QVir_Analysis/Jet Comparison.png");

    return 0;
}
