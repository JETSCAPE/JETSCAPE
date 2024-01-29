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
    int maxdir = directories.size();
    chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build");
    
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //declaring vectors to hold results
    vector<vector<double>> jetPredicts, jetPredicts2, jetPredicts3;
    vector<vector<double>> hadronPredicts, pionPredicts, kaonPredicts, protonPredicts;

    //Declaring canvas and legends
    TCanvas* cHad = new TCanvas("c1","c1",1400,1200);
    TCanvas* cPion = new TCanvas("c3","c2",1400,1200);
    TCanvas* cKaon = new TCanvas("c4","c4",1400,1200);
    TCanvas* cProton = new TCanvas("c5","c5",1400,1200);
    TCanvas* cJet = new TCanvas("c2","c2",1400,1200);
    TLegend hadronleg(.75,.35,.9,.9,"Sources");
    TLegend pionleg(.75,.35,.9,.9,"Sources");
    TLegend kaonleg(.75,.35,.9,.9,"Sources");
    TLegend protonleg(.75,.35,.9,.9,"Sources");
    TLegend jetleg(.75,.35,.9,.9,"Sources");

    //Declaring graphs
    TMultiGraph* hadrons = new TMultiGraph("Hadron Spectra","Hadron Spectra");
    TGraphErrors* jetscapehadrons[directories.size()];
    TMultiGraph* pions = new TMultiGraph("Pion Spectra","Pion Spectra");
    TGraphErrors* jetscapepions[directories.size()];
    TMultiGraph* kaons = new TMultiGraph("Kaon Spectra","Kaon Spectra");
    TGraphErrors* jetscapekaons[directories.size()];
    TMultiGraph* protons = new TMultiGraph("Proton Spectra","Proton Spectra");
    TGraphErrors* jetscapeprotons[directories.size()];
    TMultiGraph* jets = new TMultiGraph("Jet Spectra","Jet Spectra");
    TGraphErrors* jetscapejets[directories.size()];

    //hadron data file
    TFile hadron_file("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/HadronData.root");
    TDirectory* hadrondir = (TDirectory*)hadron_file.Get("Table 1");
    TGraphErrors* hadronData = (TGraphErrors*)hadrondir->Get("Graph1D_y1");
    hadronData->SetMarkerStyle(kFullDotLarge);
    hadrons->Add(hadronData,"ap");
    hadronleg.AddEntry(hadronData,"Data","ep");

    //Jet data file
    TFile jet_file("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/JetData.root");
    TDirectory* jetdir = (TDirectory*)jet_file.Get("Table 4");
    TGraphErrors* jetData = (TGraphErrors*)jetdir->Get("Graph1D_y2");
    jetData->SetMarkerStyle(kFullDotLarge);
    jets->Add(jetData,"ap");
    jetleg.AddEntry(jetData,"Data","ep");

    //ID hadron data file
    TFile idhadron_file("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/LHC-ID-hads.root");
    TDirectory* piondir = (TDirectory*)idhadron_file.Get("Table 1");
    TGraphErrors* pionData = (TGraphErrors*)piondir->Get("Graph1D_y3");
    pionData->SetMarkerStyle(kFullDotLarge);
    pions->Add(pionData,"ap");
    pionleg.AddEntry(pionData,"Data","ep");
    
    TDirectory* kaondir = (TDirectory*)idhadron_file.Get("Table 2");
    TGraphErrors* kaonData = (TGraphErrors*)kaondir->Get("Graph1D_y3");
    kaonData->SetMarkerStyle(kFullDotLarge);
    kaons->Add(kaonData,"ap");
    kaonleg.AddEntry(kaonData,"Data","ep");
    
    TDirectory* protondir = (TDirectory*)idhadron_file.Get("Table 3");
    TGraphErrors* protonData = (TGraphErrors*)protondir->Get("Graph1D_y3");
    protonData->SetMarkerStyle(kFullDotLarge);
    protons->Add(protonData,"ap");
    protonleg.AddEntry(protonData,"Data","ep");
 
    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < maxdir; i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistHadrons = (TH1D*)file.Get("hadrons"); cout << "Got hadrons. ";
        TH1D* tempHistPions = (TH1D*)file.Get("identified pions"); cout << "Got pions. ";
        TH1D* tempHistKaons = (TH1D*)file.Get("identified kaons"); cout << "Got kaons. ";
        TH1D* tempHistProtons = (TH1D*)file.Get("identified protons"); cout << "Got protons. ";
        TH1D* tempHistJets = (TH1D*)file.Get("jet radius 0.3"); cout << "Got jets 1. ";
        TH1D* tempHistJets2 = (TH1D*)file.Get("jet radius 0.2"); cout << "Got jets 2. ";
        TH1D* tempHistJets3 = (TH1D*)file.Get("jet radius 0.4"); cout << "Got jets 3." << endl;

        //adding to jets graph
        tempHistJets->Scale(1e6);
        tempHistJets2->Scale(1e6);
        tempHistJets3->Scale(1e6);
        jetscapehadrons[i] = (TGraphErrors*)histToGraph(tempHistHadrons).Clone();
        jetscapepions[i] = (TGraphErrors*)histToGraph(tempHistPions).Clone();
        jetscapekaons[i] = (TGraphErrors*)histToGraph(tempHistKaons).Clone();
        jetscapeprotons[i] = (TGraphErrors*)histToGraph(tempHistProtons).Clone();
        jetscapejets[i] = (TGraphErrors*)histToGraph(tempHistJets).Clone();

        //csv file creation
        string hadronCSVname = to_string(i) + "_hadrons.csv";
        string jetCSVname = to_string(i) + "_jets.csv";
        tempHistHadrons->GetXaxis()->SetTitle("pT (GeV)");
        tempHistJets->GetXaxis()->SetTitle("pT (GeV)");
        tempHistHadrons->GetYaxis()->SetTitle("Particle Yields (Gev^-2*c^3)");
        tempHistJets->GetYaxis()->SetTitle("Jet Differenctial Cross-section (mb/GeV/c)");
        histToCSV(tempHistHadrons, hadronCSVname);
        histToCSV(tempHistJets, jetCSVname);

        //setting graph characteristics
        color++;
        string name = to_string(i);
        jetscapehadrons[i]->SetName(name.c_str());
        jetscapepions[i]->SetName(name.c_str());
        jetscapekaons[i]->SetName(name.c_str());
        jetscapeprotons[i]->SetName(name.c_str());
        jetscapejets[i]->SetName(name.c_str());
        jetscapehadrons[i]->SetLineColor(color);
        jetscapepions[i]->SetLineColor(color);
        jetscapekaons[i]->SetLineColor(color);
        jetscapeprotons[i]->SetLineColor(color);
        jetscapejets[i]->SetLineColor(color);
        hadrons->Add(jetscapehadrons[i],"lX");
        pions->Add(jetscapepions[i],"lX");
        kaons->Add(jetscapekaons[i],"lX");
        protons->Add(jetscapeprotons[i],"lX");
        jets->Add(jetscapejets[i],"lX");
        hadronleg.AddEntry(jetscapehadrons[i], name.c_str(), "l");
        pionleg.AddEntry(jetscapepions[i], name.c_str(), "l");
        kaonleg.AddEntry(jetscapekaons[i], name.c_str(), "l");
        protonleg.AddEntry(jetscapeprotons[i], name.c_str(), "l");
        jetleg.AddEntry(jetscapejets[i], name.c_str(), "l");

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

    //Drawing hadrons
    cHad->cd();
    cHad->SetLeftMargin(0.15);
    cHad->SetLogx(); cHad->SetLogy();
    hadrons->Draw();
    hadronleg.DrawClone("Same");
    cHad->Print("QVir_Analysis/Hadron Comparison.png");
    cHad->Close();

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

    //Drawing jets
    cJet->cd();
    cJet->SetLeftMargin(0.15);
    cJet->SetLogy();
    jets->Draw();
    jetleg.DrawClone("Same");
    cJet->Print("QVir_Analysis/Jet Comparison.png");
    cJet->Close();

    //exit behavior
    hadron_file.Close();
    jet_file.Close();
    return 0;
}
