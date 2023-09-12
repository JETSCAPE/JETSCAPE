// This program reads "test_out.dat" produced from the JETSCAPE and Generates spectrums for jets and charged hadron yield.
// JETSCAPE by default gives pTHatCrossSection in milibarn and particle's momentum in GeV.
// With slight change, this program also allows reading multiple pTHatBins data and produce weighted spectra. 
// Output of this program is a ROOT file which contains following plots
// 1. Number of jets vs pT of the jet. (Graph name is CountVspTJetSpectrumBinpTHatMin_pTHatMax)
// 2. Number of charged hadron vs pT of the hadron. (Graph name is CountVspTSingleHadronSpectrumBinpTHatMin_pTHatMax)
// 3. Weighted differential jet cross section vs pT of the jet. Here Weighted differential crosssection means = sigmapTHat*dN/(dpTJet*dEta*TotalEvents)
//    (Graph name is DifferentialJetCrossSectionBinpTHatMin_pTHatMax)
//
// 4. Weighted charged hadron differential yield vs pT of hadron. Here Weighted differential Yield = sigmapTHat*dN/(TotalInelasticCrossSection*2*PI*dpTHadron*dEta*TotalEvents)  
//    (Graph name is DifferentialSingleHadronYieldBinpTHatMin_pTHatMax)
//
// 5. Total differential jet cross section (i.e. summed over all pTHat bins) vs pT of the jet. (if you have multiple pTHatBins data)
//    (Graph name is TotalDifferentialJetCrossSection)
//
// 6. Total differential charged hadron yield (i.e. summed over all pTHat bins) vs pT of the hadron.
//    (Graph name is TotalDifferentialSingleHadronYield)
//
// 7. pTHat crosssection i.e. hard scattering crosssection vs pTHat bin.
//    (Graph name is HardCrossSection) 
// 
// Note, plots are saved in two ROOT format TH1D and TGraphErrors.
// For jet spectrum, we use anti-kT algorithm with jet radius, R=0.3 and |eta_jet|<2.0. Inside jet cone, we include all particle except neutrinos (CMS def).
// For charged hadron spectrum, we use |eta_hadron|<1.0. Only charged hadrons are included in the spectrum.

// Authorship: written by Amit Kumar, revised by Shanshan Cao.

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

//my library
#include "analysis.h"

using namespace std;
using namespace Jetscape;

//using namespace Pythia8;
int main(int argc, char* argv[]){
    //change to data set directory
    string analysisDir = argv[1];
    chdir(analysisDir.c_str());
    int StartTime = time(NULL);

    // calling ROOT and pythia
    TApplication theApp("hist", &argc, argv);
    Pythia8::Pythia pythia;//("",false);
    
    //System variables
    double Ecm = 91.2;
    
    //Variables for single hadron spectrum and identified hadrons
    double SingleHadronpTBin[] = {0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.43,0.46,0.49,0.52,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1};
    int NpTSingleHadronBin = sizeof(SingleHadronpTBin)/sizeof(SingleHadronpTBin[0])-1;

    //reco hadrons
    TH1D *HistRecoHadron = new TH1D("Frag Hadrons", "Frag Hadrons", NpTSingleHadronBin, SingleHadronpTBin);
    
    //particle info declaration
    int  SN=0,PID=0;
    double Px, Py, Pz, E, Eta, Phi,pStat;
    int Events =0;
    vector<shared_ptr<Hadron>> hadrons;
    vector<double> pts, ids, recos;
    
    // Create a file on which histogram(s) can be saved.
    char outFileName[1000];
    sprintf(outFileName,"SpectraBin.root");
    TFile* outFile = new TFile( outFileName, "RECREATE");
    
    //opening run file
    char HadronFile[300], pTBinString[100];
    sprintf(HadronFile,"run.dat");
    auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
       
    //event loop
    while (!myfile->Finished()){
        myfile->Next();
        hadrons = myfile->GetHadrons();
        Events++;

        //hadron loop
        for(unsigned int i=0; i<hadrons.size(); i++){
            SN = i;
            PID= hadrons[i].get()->pid();
            E  = hadrons[i].get()->e();
            Px = hadrons[i].get()->px();
            Py = hadrons[i].get()->py();
            Pz = hadrons[i].get()->pz();
            Eta = hadrons[i].get()->eta();
            Phi = hadrons[i].get()->phi();
            pStat = hadrons[i].get()->pstat();
            double PT = TMath::Sqrt((Px*Px) + (Py*Py));
            double P = TMath::Sqrt((Px*Px) + (Py*Py) + (Pz*Pz));

            //Add to multiplicity and prepare for xp filling
            if(PT > 0.2 && abs(Eta) < 1.74 && (fabs(PID) > 100 || abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15) && pythia.particleData.charge(PID)!= 0){
                //multiplicity++;
                pts.push_back(P);
                ids.push_back(PID);
            }
        

        
        }
    }
    myfile->Close();

    //Write unprocessed data into a root file
    HistTempSingleHadron->Write();
    
    //scaling histograms
    HistMultiplicity->Scale(1/(double)Events); //times 2 to match the data
    HistTempSingleHadron->Scale(1/(double)Events,"width");
    HistRecoHadron->Scale(1/(double)Events,"width");
    
    TVector EventInfo(3);
    EventInfo[0] = 1;
    EventInfo[1] = 0;
    EventInfo[2] = Events;
    EventInfo.Write("EventInfo");
    outFile->Close();

    //ee Jets debugging
    //cout << integral(HistTempJet) << endl;

    //create root file for total plots
    TFile* totalroot = new TFile("totals.root", "RECREATE");
    HistTempSingleHadron->Write("hadrons");
    HistRecoHadron->Write("Reco Hadron Ratio");

    //ratio plots for comparison
    TGraphErrors* thrustgraph = (TGraphErrors*) thrustdir->Get("Graph1D_y1");

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;

    return 0;
}
