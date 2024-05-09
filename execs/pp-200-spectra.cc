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
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFMidPointPlugin.hh"
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

    //get list of cross sections
    double xsectotal = 42.7; //experimental value: https://arxiv.org/pdf/2005.00776.pdf
    //for(int k = 0; k<NpTHardBin; k++) cout << pTHatMin[k] << " " << pTHatMax[k] << " " << xsecList[k]*100000 << endl; //debugging line
    //cout << xsectotal << endl;

    //Cut variables
    double DetectorEtaCut= 2.6;
    double SingleHadronEtaCut = 0.5;
    double idHadronEtaCut = 0.6;
    double hardPionEtaCut = 0.35;
    double softend = 6.0;

    double jetYmin = 0.2;
    double jetYmax = 0.8;
    double jetR = 0.4;
    double overlap_threshold = 0.5;
    fastjet::CDFMidPointPlugin cdfcone(jetR, overlap_threshold);
    fastjet::JetDefinition jetDef(& cdfcone);

    //Variables for single pion spectrum
    TFile idhadron_file("/data/rjfgroup/rjf01/cameron.parker/data/PHENIX-ID-hads.root");
    TDirectory* piondir = (TDirectory*)idhadron_file.Get("Table 1");
    TH1D* piondata = (TH1D*) piondir->Get("Hist1D_y1");
    TGraphErrors* piongraph = (TGraphErrors*) piondir->Get("Graph1D_y1");
    TH1D *HistTotalPions = getBlankCopy(piondata,"Pion Spectrum","Pion Spectrum");

    //kaons
    TDirectory* kaondir = (TDirectory*)idhadron_file.Get("Table 2");
    TH1D* kaondata = (TH1D*) kaondir->Get("Hist1D_y1");
    TGraphErrors* kaongraph = (TGraphErrors*) kaondir->Get("Graph1D_y1");
    TH1D *HistTotalKaons = getBlankCopy(kaondata,"Kaon Spectrum","Kaon Spectrum");

    //protons
    TDirectory* protondir = (TDirectory*)idhadron_file.Get("Table 4");
    TH1D* protondata = (TH1D*) protondir->Get("Hist1D_y1");
    TGraphErrors* protongraph = (TGraphErrors*) protondir->Get("Graph1D_y1");
    TH1D *HistTotalProtons = getBlankCopy(protondata,"Proton Spectrum","Proton Spectrum");

    //hard pions
    double hardPionBins[] = {8.0,8.5,9.0,9.5,10.0,12.0,14.0,16.0,18.0,20.0};
    int NhardPionBins = sizeof(hardPionBins)/sizeof(hardPionBins[0])-1;
    TH1D *HistTotalHardPions = new TH1D("Hadron Spectrum", "Hadron Spectrum pT", NhardPionBins, hardPionBins);

    //Variables for single hadron spectrum
    double hadpTBin[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.0,3.35,3.8,4.4,5.1,6,7,8,9,10};
    int NpThadBin = sizeof(hadpTBin)/sizeof(hadpTBin[0])-1;
    TH1D *HistTotalHads = new TH1D("Pi0 Spectrum", "Pi0 Spectrum pT", NpThadBin, hadpTBin);

    //Jet declarations
    TFile jetFile("/data/rjfgroup/rjf01/cameron.parker/data/STAR-jets.root");
    //TDirectory* jetDir = (TDirectory*)jetFile.Get("Figure 2 - Combined HT");
    TH1D* jetdata = (TH1D*) jetFile.Get("hist");
    TGraphErrors* jetgraph = (TGraphErrors*) jetFile.Get("jets");
    TH1D* HistTotalJets = getBlankCopy(jetdata, "Total Jet Spectra", "Jet Spectra");

    //graph declaration for adding hadron spectra
    TMultiGraph* hadronComp = new TMultiGraph();
    TGraph* hadronComponents[NpTHardBin];
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat.gz", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"SpectraBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        binfiles[k] = totalroot->mkdir(outFileName);
        binfiles[k]->cd();
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int Events =0;
        int TriggeredJetNumber=0;

        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *tempPions = getBlankCopy(piondata,"Pion Spectrum", "Pion Spectrum pT");
        TH1D *tempKaons = getBlankCopy(kaondata,"Kaon Spectrum", "Kaon Spectrum pT");
        TH1D *tempProtons = getBlankCopy(protondata,"Proton Spectrum", "Proton Spectrum pT");
        TH1D *tempHads = new TH1D("Hadron Spectrum", "Hadron Spectrum pT", NpThadBin, hadpTBin);
        TH1D *tempHardPions = new TH1D("Pi0 Spectrum", "Pi0 Spectrum pT", NhardPionBins, hardPionBins);
        TH1D *tempJets = getBlankCopy(jetdata, "Bin Jet Spectra", "Jet Spectra");

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        
        //actually reading in
        while (!myfile->Finished()){
            cout << HadronFile << ": ";
            myfile->Next();
            hadrons = myfile->GetHadrons();
            vector<fastjet::PseudoJet> fjInputs;

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
                Y = hadrons[i].get()->rapidity();
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));      


                //cutting for specific regimes
                double strength = 1.0; //smoothing between smooth and hard transition 
                if(k == 0 && PT > softend)
                    strength = 0;
                if(k != 0 && PT < softend)
                    strength = 0;         

                if(fabs(Y) < idHadronEtaCut){
                    if(PID == 211) tempPions->Fill(PT, strength);
                    if(PID == 321) tempKaons->Fill(PT, strength);
                    if(PID == 2212) tempProtons->Fill(PT, strength);
                } 

                if(fabs(Eta) < hardPionEtaCut && abs(PID) == 211)
                    tempHardPions->Fill(PT, strength/2.0);
                
                // Add this particle into SingleHadron spectrum
                if(fabs(Eta) < SingleHadronEtaCut && fabs(PID)>100 &&  pythia.particleData.charge(PID)!=0){
                    tempHads->Fill(PT, strength/2.0);
                }

                //jet input prep
                if(fabs(Eta) < DetectorEtaCut && PT>0.2 && PID!=12 && PID!=14 && PID!=16 && PID!=18 )
                    fjInputs.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
            }

            //jet running
            if(k==0) continue;
            fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
            vector<fastjet::PseudoJet> UnSortedJets = clustSeq.inclusive_jets(5.0);
            for(int ijet = 0; ijet < UnSortedJets.size(); ijet++){
                if(abs(UnSortedJets[ijet].rapidity()) < jetYmax && abs(UnSortedJets[ijet].rapidity()) > jetYmin)
                    tempJets->Fill(UnSortedJets[ijet].pt(),1.0/2.0);
            }

        }

        //xsec stuff
        double HardCrossSection = myfile->GetSigmaGen();
        double HardCrossSectionError = myfile->GetSigmaErr();
        if(k == 0) xsectotal = HardCrossSection; //set for first bin to match experimental value; end of reading cross section
        
        //event count handling
        eventCount.push_back(Events);
        
        //Write histogram into a root file
        tempPions->Write();
        tempKaons->Write();
        tempProtons->Write();
        tempHads->Write();
        tempHardPions->Write();
        tempJets->Write();
        
        //add to totals histograms
        HistTotalPions->Add(tempPions,HardCrossSection/Events);
        HistTotalKaons->Add(tempKaons,HardCrossSection/Events);
        HistTotalProtons->Add(tempProtons,HardCrossSection/Events);
        HistTotalHads->Add(tempHads,HardCrossSection/(xsectotal*Events));
        HistTotalHardPions->Add(tempHardPions,HardCrossSection/Events);
        HistTotalJets->Add(tempJets,HardCrossSection/Events);
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");
        totalroot->cd();
    } //k-loop ends here (pTHatBin loop)

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    scaleBins(HistTotalPions,(1.0/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalKaons,(1.0/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalProtons,(1.0/(2*M_PI*2.0*idHadronEtaCut)));
    scaleBins(HistTotalHads,(1.0/(2*M_PI*2.0*SingleHadronEtaCut)));
    scaleBins(HistTotalHardPions,(1.0/(2*M_PI*2.0*hardPionEtaCut)));
    HistTotalJets->Scale((1.0e9/(2*M_PI*(jetYmax-jetYmin))),"width");

    //Plotting
    myRatioPlot(piongraph, HistTotalPions, "Pion Yields", true, true);
    myRatioPlot(kaongraph, HistTotalKaons, "Kaon Yields", true, true);
    myRatioPlot(protongraph, HistTotalProtons, "Proton Yields", true, true);
    myRatioPlot(jetgraph, HistTotalJets, "Jet Yields", false, true);
 	
    //create root file for total plots
    HistTotalPions->Write("raw pions"); smoothBins(HistTotalPions); HistTotalPions->Write("identified pions");
    HistTotalKaons->Write("raw kaons"); smoothBins(HistTotalKaons); HistTotalKaons->Write("identified kaons");
    HistTotalProtons->Write("raw protons"); smoothBins(HistTotalProtons); HistTotalProtons->Write("identified protons");
    HistTotalHads->Write("raw hads"); smoothBins(HistTotalHads); HistTotalHads->Write("identified hads");
    HistTotalHardPions->Write("raw hard pions"); smoothBins(HistTotalHardPions); HistTotalHardPions->Write("identified hard pions");
    HistTotalJets->Write("jets");
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
