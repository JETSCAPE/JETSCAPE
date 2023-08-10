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
    chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build");

    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //array to hold values
    vector<vector<double>> pionPredicts;
    vector<vector<double>> hadPredicts;

    //Declaring graphs
    TCanvas* cHad = new TCanvas("c1","c1",1400,1200);
    TLegend hadronleg(.75,.35,.9,.9,"Sources");
    TMultiGraph* hadrons = new TMultiGraph("Hadron Spectra","Hadron Spectra");
    TGraphErrors* jetscapehadrons[directories.size()];

    //hadron data file
    TFile hadron_file("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/STARhadData.root");
    TDirectory* hadrondir = (TDirectory*)hadron_file.Get("Table 1");
    TGraphErrors* hadronData = (TGraphErrors*)hadrondir->Get("Graph1D_y8");
    hadronData->SetMarkerStyle(kFullDotLarge);
    hadrons->Add(hadronData,"ap");
    hadronleg.AddEntry(hadronData,"Data","ep");
 
    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < directories.size(); i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistPions = (TH1D*)file.Get("identified pions"); cout << "Got pions. ";
        TH1D* tempHistHads = (TH1D*)file.Get("identified hads"); cout << "Got hadrons. " << endl;

        //dat file data for Bayes analysis
        vector<double> piontemp, hadtemp;
        for(int j = 1; j <= tempHistPions->GetNbinsX(); j++) piontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistHads->GetNbinsX(); j++) hadtemp.push_back(tempHistHads->GetBinContent(j));
        pionPredicts.push_back(piontemp);
        hadPredicts.push_back(hadtemp);
        jetscapehadrons[i] = (TGraphErrors*)histToGraph(tempHistHads).Clone();

        //setting graph characteristics
        color++;
        string name = to_string(i);
        jetscapehadrons[i]->SetName(name.c_str());
        jetscapehadrons[i]->SetLineColor(color);
        hadrons->Add(jetscapehadrons[i],"lX");
        hadronleg.AddEntry(jetscapehadrons[i], name.c_str(), "l");
        
        file.Close();
    }

    //making said dat files
    chdir(input.c_str());
    char tmp[256]; getcwd(tmp, 256); cout << "Current working directory: " << tmp << endl;
    makeDatFile(pionPredicts, "PionSpectraPrediction", "# Version 1.0\n# Pion Spectra for parameters.txt");
    makeDatFile(hadPredicts, "HadronSpectraPrediction", "# Version 1.0\n# Hadron Spectra for parameters.txt");
    
    //Drawing hadrons
    cHad->cd();
    cHad->SetLeftMargin(0.15);
    cHad->SetLogx(); cHad->SetLogy();
    hadrons->Draw();
    hadronleg.DrawClone("Same");
    cHad->Print("QVir_Analysis/Hadron Comparison.png");
    cHad->Close();

    //exit behavior
    hadron_file.Close();
    return 0;
}
