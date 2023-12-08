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
#include "THStack.h"

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
    
    //Variables for single hadron spectrum and identified hadrons from https://arxiv.org/pdf/hep-ex/9909032.pdf
    double SingleHadronpTBin[] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
    double spectrumData[] = {7.47,9.03,10.42,10.76,9.89,8.97,8.17,6.94,6.73,5.56,4.94,3.49,3.13,2.00,1.27,0.50,0.27,0.06};
    double spectrumError[] = {0.63,0.49,0.44,0.43,0.38,0.35,0.32,0.28,0.27,0.24,0.22,0.18,0.17,0.14,0.11,0.07,0.05,0.03};
    int NpTSingleHadronBin = sizeof(SingleHadronpTBin)/sizeof(SingleHadronpTBin[0])-1;
    double spectrumX[NpTSingleHadronBin], spectrumWidth[NpTSingleHadronBin];
    for(int i=0; i<NpTSingleHadronBin; i++){
        spectrumX[i] = (SingleHadronpTBin[i+1]+SingleHadronpTBin[i])/2;
        spectrumWidth[i] = (SingleHadronpTBin[i+1]-SingleHadronpTBin[i])/2;
    }

    //reco hadrons
    TH1D *HistDstars = new TH1D("Frag Hadrons", "Frag Hadrons", NpTSingleHadronBin, SingleHadronpTBin);
    TGraphErrors *data = new TGraphErrors(NpTSingleHadronBin, spectrumX, spectrumData, spectrumWidth, spectrumError);
    
    //particle info declaration
    int  SN=0,PID=0;
    double Px, Py, Pz, E, Eta, Phi,pStat;
    int Events=0;
    vector<shared_ptr<Hadron>> hadrons;
    vector<double> pts, ids, recos;
    
    // Create a file on which histogram(s) can be saved.
    char outFileName[1000];
    sprintf(outFileName,"SpectraBin.root");
    TFile* outFile = new TFile(outFileName, "RECREATE");
    
    //opening run file
    char HadronFile[300], pTBinString[100];
    sprintf(HadronFile,"run.dat.gz");
    auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
       
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
            if(abs(PID) == 413){
                //multiplicity++;
                double xe = 2*E/Ecm;
                HistDstars->Fill(xe);
            }
        }
    }
    myfile->Close();
    
    //scaling histograms
    HistDstars->Scale((1000.0*0.0389*0.677)/(double)Events,"width");
    HistDstars->Write("D stars");
    ratioPlot(data,HistDstars,"Dstarxe","x_{e}","1/N_{events} dN_{D*}/d_{xe}",false,false);
    
    TVector EventInfo(3);
    EventInfo[0] = 1;
    EventInfo[1] = 0;
    EventInfo[2] = Events;
    EventInfo.Write("EventInfo");
    outFile->Close();

    HistDstars->Print("all");
    data->Print("all");

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;

    return 0;
}
