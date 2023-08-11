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
#include "TProfile.h"

#include "analysis.cc"

using namespace std;
using namespace Jetscape;

//using namespace Pythia8;
int main(int argc, char* argv[]){
    //change to data set directory
    string analysisDir = argv[1];
    chdir(analysisDir.c_str());

    int nListJets =1;
    int StartTime = time(NULL);
    // Create the ROOT application environment and pythia.
    TApplication theApp("hist", &argc, argv);
    Pythia8::Pythia pythia;//("",false);
    
    //Total analysis variables
    //Reading ptHat bins from list of made directories
    vector<vector<string>> tempvec = getDatBounds("./dat");
    vector<string> pTHatMin = tempvec[0];
    vector<string> pTHatMax = tempvec[1];
    int NpTHardBin = pTHatMin.size();
    //for(int i = 0; i < pTHatMin.size(); i++) cout << pTHatMin[i] << endl; //debugging line
    vector<int> eventCount;

    //xsec total running count
    double xsectotal = 0.0;

    //Cut variables
    double idHadronEtaCut = 1;
    double softend = 6.0;
    
    //Variables for single pion spectrum
    int NpTpionBin = 22;
    double pionpTBin[NpTpionBin]; for(int i=0; i<=NpTpionBin; i++) pionpTBin[i] = (0.05*i)+0.1;
    TH1D *HistTotalPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin); //identified hadrons hists
    
    //Variables for single kaon spectrum
    int NpTkaonBin = 17;
    double kaonpTBin[NpTkaonBin]; for(int i=0; i<=NpTkaonBin; i++) kaonpTBin[i] = (0.05*i)+0.2;
    TH1D *HistTotalKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", NpTkaonBin, kaonpTBin); //identified hadrons hists
    
    //Variables for single proton spectrum
    int NpTprotonBin = 27;
    double protonpTBin[NpTprotonBin]; for(int i=0; i<=NpTprotonBin; i++) protonpTBin[i] = (0.05*i)+0.35;
    TH1D *HistTotalProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protonpTBin); //identified hadrons hists
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass;
        int Events =0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s.root",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        TFile* outFile = new TFile( outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *tempPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);
        TH1D *tempKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", NpTkaonBin, kaonpTBin);
        TH1D *tempProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protonpTBin);

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        
        //actually reading in
        while (!myfile->Finished()){
            myfile->Next();
            hadrons = myfile->GetHadrons();

            //cout<<"Number of hadrons is: " << hadrons.size() << endl;
            Events++;
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
                mass = hadrons[i].get()->restmass();
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));      

                //cutting for specific regimes
                if(k == 0 && PT > softend)
                    continue;
                if(k != 0 && PT < softend)
                    continue;

                double strength = 1; //smoothing between smooth and hard transition          

                if(fabs(Eta) < idHadronEtaCut){
                    if(abs(PID) == 211) tempPions->Fill(PT, strength/2);
                    if(abs(PID) == 321) tempKaons->Fill(PT, strength/2);
                    if(abs(PID) == 2212) tempProtons->Fill(PT, strength/2);
                } 
            }
        }

        //xsec stuff
        double HardCrossSection = myfile->GetSigmaGen();
        double HardCrossSectionError =  myfile->GetSigmaErr();
        xsectotal += HardCrossSection; //set for first bin to match experimental value; end of reading cross section
        
        //event count handling
        eventCount.push_back(Events);
        
        //Write histogram into a root file
        tempPions->Write();
        tempKaons->Write();
        tempProtons->Write();
        
        //add to totals histograms
        HistTotalPions->Add(tempPions,HardCrossSection/Events);
        HistTotalKaons->Add(tempKaons,HardCrossSection/Events);
        HistTotalProtons->Add(tempProtons,HardCrossSection/Events);
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");

        outFile->Close();
    } //k-loop ends here (pTHatBin loop)

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    scaleBins(HistTotalPions,(1/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalKaons,(1/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalProtons,(1/(2*M_PI*2.0*idHadronEtaCut)));
 	
    //create root file for total plots
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");
    HistTotalPions->Write("raw pions"); smoothBins(HistTotalPions); HistTotalPions->Write("identified pions");
    HistTotalPions->Write("raw kaons"); smoothBins(HistTotalKaons); HistTotalPions->Write("identified kaons");
    HistTotalPions->Write("raw protons"); smoothBins(HistTotalProtons); HistTotalPions->Write("identified protons");
    totalroot->Close();

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    
    //test comment
    //debugging
    //for(int i = 0; i < NpTJetBin; i++) cout << DifferentialJetTotal[i]*1000000 << " " << DifferentialJetTotalErrors[i]*1000000 << endl;
    //for(int i = 0; i < NpTSingleHadronBin; i++) cout << DifferentialHadronTotal[i] << " " << DifferentialHadronTotalErrors[i]*1000000 << endl;
    //for(int i = 0; i < NpTHardBin; i++) cout << pTHardBinSingleHadronBinError << endl;
    //for(int i = 0; i < eventCount.size(); i++) cout << pTHatMin[i] << " " << pTHatMax[i] << " " << eventCount[i] << " " << xsecList[i] << endl;
    return 0;
}
