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
#include "TDirectory.h"

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
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");
    Pythia8::Pythia pythia;//("",false);
    
    //Total analysis variables
    //Reading ptHat bins from list of made directories
    vector<vector<string>> tempvec = getDatBounds("./dat");
    vector<string> pTHatMin = tempvec[0];
    vector<string> pTHatMax = tempvec[1];
    int NpTHardBin = pTHatMin.size();
    //for(int i = 0; i < pTHatMin.size(); i++) cout << pTHatMin[i] << endl; //debugging line
    vector<int> eventCount;
    TDirectory* binfiles[NpTHardBin];

    //xsec total running count
    double xsectotal = 0.0;

    //Cut variables
    double idHadronYCut = 0.5;
    double softend = 6.0;
    
    //reading data to get bins
    TFile dataroot( "/data/rjfgroup/rjf01/cameron.parker/data/LHC13000.root");
    TDirectory* piondir = (TDirectory*)dataroot.Get("Table 1"); TH1D* piondata = (TH1D*)piondir->Get("Hist1D_y1");
    TDirectory* kaondir = (TDirectory*)dataroot.Get("Table 2"); TH1D* kaondata = (TH1D*)kaondir->Get("Hist1D_y1");
    TDirectory* protondir = (TDirectory*)dataroot.Get("Table 6"); TH1D* protondata = (TH1D*)protondir->Get("Hist1D_y1");

    //Variables for ID hadron hists
    TH1D *HistTotalPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat.gz", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Y, Phi, pStat, mass;
        int Events =0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"SpectraBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        binfiles[k] = totalroot->mkdir(outFileName);
        binfiles[k]->cd();
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *tempPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        
        //actually reading in
        while (!myfile->Finished()){
            cout << HadronFile << ": ";
            try{
                myfile->Next();
                hadrons = myfile->GetHadrons();
            }
            catch(...){
                break;
            }

            //cout<<"Number of hadrons is: " << hadrons.size() << endl;
            Events++;
            for(unsigned int i=0; i<hadrons.size(); i++){
                SN = i;
                PID= hadrons[i].get()->pid();
                E  = hadrons[i].get()->e();
                Px = hadrons[i].get()->px();
                Py = hadrons[i].get()->py();
                Pz = hadrons[i].get()->pz();
                Y = hadrons[i].get()->rapidity();
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

                if(fabs(Y) < idHadronYCut){
                    if(abs(PID) == 211) tempPions->Fill(PT, strength/2.);
                    if(abs(PID) == 321) tempKaons->Fill(PT, strength/2.);
                    if(abs(PID) == 2212) tempProtons->Fill(PT, strength/2.);
                } 
            }
        }

        //xsec stuff
        double HardCrossSection = myfile->GetSigmaGen();
        double HardCrossSectionError =  myfile->GetSigmaErr();
        if(k == 0) xsectotal = HardCrossSection; //set for first bin to match experimental value; end of reading cross section
        
        //event count handling
        eventCount.push_back(Events);
        
        //Write histogram into a root file
        tempPions->Write();
        tempKaons->Write();
        tempProtons->Write();
        
        //add to totals histograms
        HistTotalPions->Add(tempPions,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalKaons->Add(tempKaons,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalProtons->Add(tempProtons,HardCrossSection/(1.0*Events*xsectotal));
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");

        totalroot->cd();
    } //k-loop ends here (pTHatBin loop)

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    HistTotalPions->Scale(1/(2.0*idHadronYCut),"width");
    HistTotalKaons->Scale(1/(2.0*idHadronYCut),"width");
    HistTotalProtons->Scale(1/(2.0*idHadronYCut),"width");
 	
    //create root file for total plots
    HistTotalPions->Write("raw pions"); smoothBins(HistTotalPions); HistTotalPions->Write("identified pions");
    HistTotalKaons->Write("raw kaons"); smoothBins(HistTotalKaons); HistTotalKaons->Write("identified kaons");
    HistTotalProtons->Write("raw protons"); smoothBins(HistTotalProtons); HistTotalProtons->Write("identified protons");
    totalroot->Close();
    HistTotalPions->Print("all");

    //hadron graphs
    myRatioPlot((TGraphErrors*)piondir->Get("Graph1D_y1"), HistTotalPions, "Pion Yields", true, true);
    myRatioPlot((TGraphErrors*)kaondir->Get("Graph1D_y1"), HistTotalKaons, "Kaon Yields", true, true);
    myRatioPlot((TGraphErrors*)protondir->Get("Graph1D_y1"), HistTotalProtons, "Proton Yields", true, true);

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
