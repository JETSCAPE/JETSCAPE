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
    int maxdir = directories.size();
    chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build");

    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);

    //array to hold values
    vector<vector<double>> pionPredicts;
    vector<vector<double>> kaonPredicts;
    vector<vector<double>> protonPredicts;
    vector<vector<double>> hardPionPredicts;
    vector<vector<double>> jetPredicts;
 
    //grabbing graphs from each jetscape file
    int color = 1;
    for(int i = 0; i < maxdir; i++){
        //reading graph
        string filename = directories[i] + "/root/totals.root";
        cout << filename << ": ";
        TFile file(filename.c_str());
        TH1D* tempHistPions = (TH1D*)file.Get("identified pions"); cout << "Got pions. ";
        TH1D* tempHistKaons = (TH1D*)file.Get("identified kaons"); cout << "Got kaons. ";
        TH1D* tempHistProtons = (TH1D*)file.Get("identified protons"); cout << "Got potons. ";
        TH1D* tempHistHardPions = (TH1D*)file.Get("hard pions"); cout << "Got hard pions. ";
        TH1D* tempHistJets = (TH1D*)file.Get("jets"); cout << "Got jets. " << endl;

        //dat file data for Bayes analysis
        vector<double> piontemp, kaontemp, protontemp, hardpiontemp, jettemp;
        for(int j = 1; j <= tempHistPions->GetNbinsX(); j++) piontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistKaons->GetNbinsX(); j++) kaontemp.push_back(tempHistKaons->GetBinContent(j));
        for(int j = 1; j <= tempHistProtons->GetNbinsX(); j++) protontemp.push_back(tempHistProtons->GetBinContent(j));
        for(int j = 1; j <= tempHistHardPions->GetNbinsX(); j++) hardpiontemp.push_back(tempHistPions->GetBinContent(j));
        for(int j = 1; j <= tempHistJets->GetNbinsX(); j++) jettemp.push_back(tempHistJets->GetBinContent(j));
        pionPredicts.push_back(piontemp);
        kaonPredicts.push_back(kaontemp);
        protonPredicts.push_back(protontemp);
        hardPionPredicts.push_back(hardpiontemp);
        jetPredicts.push_back(jettemp);
        
        file.Close();
    }

    //making said dat files
    chdir(input.c_str());
    char tmp[256]; getcwd(tmp, 256); cout << "Current working directory: " << tmp << endl;
    makeDatFile(pionPredicts, "PionSpectraPrediction", "# Version 1.0\n# Pion Spectra for parameters.txt");
    makeDatFile(kaonPredicts, "KaonSpectraPrediction", "# Version 1.0\n# Kaon Spectra for parameters.txt");
    makeDatFile(protonPredicts, "ProtonSpectraPrediction", "# Version 1.0\n# Proton Spectra for parameters.txt");
    makeDatFile(hardPionPredicts, "HardPionSpectraPrediction", "# Version 1.0\n# Hard Pion Spectra for parameters.txt");
    makeDatFile(jetPredicts, "JetSpectraPrediction", "# Version 1.0\n# Jet Spectra for parameters.txt");
    
    //exit behavior
    return 0;
}
