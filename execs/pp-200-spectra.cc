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
    if(stoi(getEcm()) != 200){
        cout << "Center of mass energy does not match." << endl;
        return 0;
    }

    //get list of cross sections
    vector<vector<double>> xsecout = getAllXsecs(pTHatMin,pTHatMax);
	vector<double> xsecList = xsecout[0];
    vector<double> xsecErrorList = xsecout[0];
    double xsectotal = 42.7; //experimental value: https://arxiv.org/pdf/2005.00776.pdf
    //for(int k = 0; k<NpTHardBin; k++) cout << pTHatMin[k] << " " << pTHatMax[k] << " " << xsecList[k]*100000 << endl; //debugging line
    //cout << xsectotal << endl;

    //Cut variables
    double DetectorEtaCut= 2.6;
    double SingleHadronEtaCut = 0.5;
    double idHadronEtaCut = 0.35;
    double softend = 6.0;
    
    //Variables for single pion spectrum
    double pionpTBin[] = {0.5,0.7,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,12,14,16,18,20};
    int NpTpionBin = sizeof(pionpTBin)/sizeof(pionpTBin[0])-1;
    TH1D *HistTotalPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin); //identified hadrons hists
    TProfile *ProfileTotalPions = new TProfile("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);

    //Variables for single hadron spectrum
    double hadpTBin[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.0,3.35,3.8,4.4,5.1,6,7,8,9,10};
    int NpThadBin = sizeof(hadpTBin)/sizeof(hadpTBin[0])-1;
    TH1D *HistTotalHads = new TH1D("Hadron Spectrum", "Hadron Spectrum pT", NpThadBin, hadpTBin);
    TProfile *ProfileTotalHads = new TProfile("Hadron Spectrum", "Hadron Spectrum pT", NpThadBin, hadpTBin);

    //graph declaration for adding hadron spectra
    TMultiGraph* hadronComp = new TMultiGraph();
    TGraph* hadronComponents[NpTHardBin];
    
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
        int TriggeredJetNumber=0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s.root",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        TFile* outFile = new TFile( outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *tempPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);
        TH1D *tempHads = new TH1D("Hadron Spectrum", "Hadron Spectrum pT", NpThadBin, hadpTBin);
        TProfile *tempPionsProf = new TProfile("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);
        TProfile *tempHadsProf = new TProfile("Hadron Spectrum", "Hadron Spectrum pT", NpThadBin, hadpTBin);

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;

        //xsec stuff
        double HardCrossSection = xsecList[k];
        double HardCrossSectionError = xsecErrorList[k];
        if(k == 0) HardCrossSection = xsectotal; //set for first bin to match experimental value; end of reading cross section
        
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

                if(fabs(Eta) < idHadronEtaCut && fabs(Phi) < pi){
                    if(abs(PID) == 211) tempPions->Fill(PT, strength/2);
                    if(abs(PID) == 211) tempPionsProf->Fill(PT, PT);
                } 
                
                // Add this particle into SingleHadron spectrum
                if(fabs(Eta) < SingleHadronEtaCut && fabs(PID)>100 &&  pythia.particleData.charge(PID)!=0){
                    tempHads->Fill(PT, strength/2);
                    tempHadsProf->Fill(PT, PT);
                }
            }
        }
        
        //event count handling
        eventCount.push_back(Events);
        
        //Write histogram into a root file
        tempPions->Write();
        tempHads->Write();
        tempPionsProf->Write();
        tempHadsProf->Write();
        
        //add to totals histograms
        HistTotalPions->Add(tempPions,HardCrossSection/Events);
        HistTotalHads->Add(tempHads,HardCrossSection/(xsectotal*Events));
        ProfileTotalPions->Add(tempPionsProf,HardCrossSection/Events);
        ProfileTotalHads->Add(tempHadsProf,HardCrossSection/(xsectotal*Events));
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");

        outFile->Close();
    } //k-loop ends here (pTHatBin loop)

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    scaleBins(HistTotalPions,ProfileTotalPions,(1/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalHads,ProfileTotalHads,(1/(2*M_PI*2.0*SingleHadronEtaCut)));
 	
    //create root file for total plots
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");
    HistTotalPions->Write("raw pions"); smoothBins(HistTotalPions); HistTotalPions->Write("identified pions");
    HistTotalHads->Write("raw hads"); smoothBins(HistTotalHads); HistTotalHads->Write("identified hads");
    ProfileTotalPions->Write("pions profile");
    ProfileTotalHads->Write("hadrons profile");
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
