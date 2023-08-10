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
    string input = argv[1]; 
    
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //declaring vectors to hold results
    vector<vector<double>> pionPredicts, kaonPredicts, protonPredicts;

    //Declaring canvas and legends
    TCanvas* cPion = new TCanvas("c3","c2",1400,1200);
    TCanvas* cKaon = new TCanvas("c4","c4",1400,1200);
    TCanvas* cProton = new TCanvas("c5","c5",1400,1200);
    TLegend pionleg(.75,.35,.9,.9,"Sources");
    TLegend kaonleg(.75,.35,.9,.9,"Sources");
    TLegend protonleg(.75,.35,.9,.9,"Sources");

    //Declaring graphs
    TMultiGraph* pions = new TMultiGraph("Pion Spectra","Pion Spectra");
    TGraphErrors* jetscapepions[directories.size()];
    TMultiGraph* kaons = new TMultiGraph("Kaon Spectra","Kaon Spectra");
    TGraphErrors* jetscapekaons[directories.size()];
    TMultiGraph* protons = new TMultiGraph("Proton Spectra","Proton Spectra");
    TGraphErrors* jetscapeprotons[directories.size()];

    //ID hadron data file
    TFile idhadron_file("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/LHC900.root");
    TDirectory* dir = (TDirectory*)idhadron_file.Get("Table 1");
    TGraphErrors* pionData = (TGraphErrors*)dir->Get("Graph1D_y1");
    pionData->SetMarkerStyle(kFullDotLarge);
    pions->Add(pionData,"ap");
    pionleg.AddEntry(pionData,"Data","ep");
    
    TGraphErrors* kaonData = (TGraphErrors*)dir->Get("Graph1D_y2");
    kaonData->SetMarkerStyle(kFullDotLarge);
    kaons->Add(kaonData,"ap");
    kaonleg.AddEntry(kaonData,"Data","ep");
    
    TGraphErrors* protonData = (TGraphErrors*)dir->Get("Graph1D_y3");
    protonData->SetMarkerStyle(kFullDotLarge);
    protons->Add(protonData,"ap");
    protonleg.AddEntry(protonData,"Data","ep");
 
    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < directories.size(); i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistPions = (TH1D*)file.Get("identified pions"); cout << "Got pions. ";
        TH1D* tempHistKaons = (TH1D*)file.Get("identified kaons"); cout << "Got kaons. ";
        TH1D* tempHistProtons = (TH1D*)file.Get("identified protons"); cout << "Got protons. " << endl;

        //adding to jets graph
        jetscapepions[i] = (TGraphErrors*)histToGraph(tempHistPions).Clone();
        jetscapekaons[i] = (TGraphErrors*)histToGraph(tempHistKaons).Clone();
        jetscapeprotons[i] = (TGraphErrors*)histToGraph(tempHistProtons).Clone();

        //setting graph characteristics
        color++;
        string name = to_string(i);
        jetscapepions[i]->SetName(name.c_str());
        jetscapekaons[i]->SetName(name.c_str());
        jetscapeprotons[i]->SetName(name.c_str());
        jetscapepions[i]->SetLineColor(color);
        jetscapekaons[i]->SetLineColor(color);
        jetscapeprotons[i]->SetLineColor(color);
        pions->Add(jetscapepions[i],"lX");
        kaons->Add(jetscapekaons[i],"lX");
        protons->Add(jetscapeprotons[i],"lX");
        pionleg.AddEntry(jetscapepions[i], name.c_str(), "l");
        kaonleg.AddEntry(jetscapekaons[i], name.c_str(), "l");
        protonleg.AddEntry(jetscapeprotons[i], name.c_str(), "l");

        //dat file data for Bayes analysis
        vector<double> piontemp, kaontemp, protontemp;
        for(int j = 1; j <= tempHistPions->GetNbinsX(); j++) piontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistKaons->GetNbinsX(); j++) kaontemp.push_back(tempHistKaons->GetBinContent(j));
        for(int j = 1; j <= tempHistProtons->GetNbinsX(); j++) protontemp.push_back(tempHistProtons->GetBinContent(j));
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        
        file.Close();
    }

    //making said dat files
    chdir(input.c_str());
    char tmp[256]; getcwd(tmp, 256); cout << "Current working directory: " << tmp << endl;
    makeDatFile(pionPredicts, "PionSpectraPrediction", "# Version 1.0\n# Pion Spectra for parameters.txt");
    makeDatFile(kaonPredicts, "KaonSpectraPrediction", "# Version 1.0\n# Kaon Spectra for parameters.txt");
    makeDatFile(protonPredicts, "ProtonSpectraPrediction", "# Version 1.0\n# Proton Spectra for parameters.txt");

    //Drawing hadrons
    cPion->cd();
    cPion->SetLeftMargin(0.15);
    cPion->SetLogx(); cPion->SetLogy();
    pions->Draw();
    pionleg.DrawClone("Same");
    cPion->Print("QVir_Analysis/Pion Comparison.png");
    cPion->Close();

    //Drawing hadrons
    cKaon->cd();
    cKaon->SetLeftMargin(0.15);
    cKaon->SetLogx(); cKaon->SetLogy();
    kaons->Draw();
    kaonleg.DrawClone("Same");
    cKaon->Print("QVir_Analysis/Kaon Comparison.png");
    cKaon->Close();

    //Drawing hadrons
    cProton->cd();
    cProton->SetLeftMargin(0.15);
    cProton->SetLogx(); cProton->SetLogy();
    protons->Draw();
    protonleg.DrawClone("Same");
    cProton->Print("QVir_Analysis/Proton Comparison.png");
    cProton->Close();

    //exit behavior
    idhadron_file.Close();
    return 0;
}
