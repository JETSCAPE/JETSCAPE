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
    
    //Variables for jet spectrum
    double JetpTBin[37]; //in GeV
    for(int i=10; i<=46; i++) JetpTBin[i-10]=i;
    int NpTJetBin = sizeof(JetpTBin)/sizeof(JetpTBin[0])-1;
    double JetRadius = 0.4;

    //dijet
    double diJetpTBin[27]; //in GeV
    for(int i=20; i<=46; i++) diJetpTBin[i-20]=i;
    int NpTdiJetBin = sizeof(diJetpTBin)/sizeof(diJetpTBin[0])-1;
    
    //Variables for single hadron spectrum and identified hadrons
    double SingleHadronpTBin[] = {0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.43,0.46,0.49,0.52,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1};
    int NpTSingleHadronBin = sizeof(SingleHadronpTBin)/sizeof(SingleHadronpTBin[0])-1;

    double pionpTBin[] = {0.005,0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.01,0.011,0.012,0.013,0.014,0.016,0.018,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.4,0.6,0.8};
    int NpTpionBin = sizeof(pionpTBin)/sizeof(pionpTBin[0])-1;

    double kaonpTBin[] = {0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.013,0.014,0.016,0.018,0.07,0.075,0.08,0.085,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.4,0.6,0.8};
    int NpTkaonBin = sizeof(kaonpTBin)/sizeof(kaonpTBin[0])-1;

    double protonpTBin[] = {0.01,0.011,0.012,0.013,0.014,0.016,0.018,0.024,0.026,0.028,0.07,0.075,0.08,0.085,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.4,0.6,0.8};
    int NpTprotonBin = sizeof(protonpTBin)/sizeof(protonpTBin[0])-1;

    //multiplicity variables
    int NMultBin = 26;
    double HadronMultiplicityBin[NMultBin];
    for(int i = 0; i <= NMultBin; i++) HadronMultiplicityBin[i] = 3+(2*i);

    //variables for thrust
    double ThrustBin[] = {0,.005,.01,.015,.02,.025,.03,.035,.04,.05,.06,.08,.1,.12,.14,.16,.18,.2,.25,.3,.35,.4};
    int NThrustBin = sizeof(ThrustBin)/sizeof(ThrustBin[0]) - 1;
    int thrustCount = 0;
    
    // Histograms for event variables
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempdiJet = new TH1D("diJetSpectrumBin", "diJet Spectrum pT", NpTdiJetBin, diJetpTBin);
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistMultiplicity = new TH1D("Charged Particle Multiplicity","Charged Particle Multiplicity",NMultBin,HadronMultiplicityBin);
    TH1D *HistThrust = new TH1D("Thrust", "Event Thrust", NThrustBin, ThrustBin);

    //reco hadrons
    TH1D *HistRecoHadron = new TH1D("Frag Hadrons", "Frag Hadrons", NpTSingleHadronBin, SingleHadronpTBin);

    //temp hists for identified hadrons
    TH1D *tempPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);
    TH1D *tempKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", NpTkaonBin, kaonpTBin);
    TH1D *tempProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protonpTBin);
    
    //jet stuff declaration
    std::vector <fjcore::PseudoJet> fjInputs;
    std::vector <int> chargeList;
    fjcore::JetDefinition jetDef(fjcore::ee_genkt_algorithm, JetRadius, -1);
    
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
    sprintf(HadronFile,"run.dat.gz");
    auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
       
    //event loop
    while (!myfile->Finished()){
        myfile->Next();
        hadrons = myfile->GetHadrons();
        int multiplicity = 0;
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

            //premature mult count
            if(pythia.particleData.charge(PID) != 0) multiplicity++;

            //jet inputs
            if(PID!=12 && PID!=14 && PID!=16 && PID!=18){
                fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                chargeList.push_back(pythia.particleData.charge(PID));
            }

            //reco flag
            //cout << hadrons[i].get()->plabel() << endl;
            if(hadrons[i].get()->plabel() <= 12)
                recos.push_back(true);
            else
                recos.push_back(false);

            //Add to multiplicity and prepare for xp filling
            if(PT > 0.2 && abs(Eta) < 1.74 && (fabs(PID) > 100 || abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15) && pythia.particleData.charge(PID)!= 0){
                //multiplicity++;
                pts.push_back(P);
                ids.push_back(PID);
            }
        }

        //mult filling
        HistMultiplicity->Fill(multiplicity);

        //Jets
        // Run Fastjet algorithm and sort jets in pT order.
        vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
        double totalE = 0;
        fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
        UnSortedJets = clustSeq.inclusive_jets();
        SortedJets = sorted_by_E(UnSortedJets);
        for (int i = 0; i < SortedJets.size(); ++i){
            double theta = 2*atan(exp(SortedJets[i].eta()));
            totalE += SortedJets[i].E();
            if(theta > pi/5 && theta < 4*pi/5) {
                HistTempJet->Fill(SortedJets[i].E());
                if(i<2) HistTempdiJet->Fill(SortedJets[i].E());
            }
            
            //debugging statements
            //cout << "jet: " << SortedJets[i].E() << endl;
            //vector <fjcore::PseudoJet> constituents = SortedJets[i].constituents();
            //for(int j = 0; j < constituents.size(); j++) cout << "\tparticle: " << constituents[j].E() << endl;
        }
        //cout << "Total Energy: " << totalE << endl;

        //Sphericity and thrust calculations
        vector<double> thrustsphericity = getThrustSphericity(hadrons);
        double thrust = thrustsphericity[0];
        double sphericity = thrustsphericity[1];

        //hadron event cuts
        if(thrust >= 0){
            //xP filling
            //HistMultiplicity->Fill(multiplicity);
            
            for(int k = 0; k < pts.size(); k++){
                double xp = 2*pts[k]/Ecm;
                HistTempSingleHadron->Fill(xp);
                if(recos[k]) HistRecoHadron->Fill(xp);
                if(abs(ids[k]) == 211) tempPions->Fill(xp);
                if(abs(ids[k]) == 321) tempKaons->Fill(xp);
                if(abs(ids[k]) == 2212) tempProtons->Fill(xp);
            }

            //other observable filling
            HistThrust->Fill(1-thrust);
            thrustCount++;
        }

        //clearing event vectors
        fjInputs.resize(0); 
        chargeList.resize(0);
        SortedJets.resize(0);
        UnSortedJets.resize(0);
        pts.resize(0);
        ids.resize(0);
        recos.resize(0);
    }
    myfile->Close();

    //Write unprocessed data into a root file
    HistMultiplicity->Write();
    HistTempSingleHadron->Write();
    HistRecoHadron->Write();
    tempPions->Write();
    tempKaons->Write();
    tempProtons->Write();
    HistThrust->Write();
    HistTempJet->Write();
    HistTempdiJet->Write();
    
    //scaling histograms
    HistMultiplicity->Scale(1.0/(double)Events); //times 2 to match the data
    HistTempSingleHadron->Scale(1.0/(double)Events,"width");
    HistRecoHadron->Scale(1.0/(double)Events,"width");
    tempPions->Scale(1.0/(double)Events,"width");
    tempKaons->Scale(1.0/(double)Events,"width");
    tempProtons->Scale(1.0/(double)Events,"width");
    HistThrust->Scale(1.0/(double)thrustCount,"width");
    HistTempJet->Scale(1.0/(double)Events);
    HistTempdiJet->Scale(1.0/(double)Events);
    
    TVector EventInfo(3);
    EventInfo[0] = 1;
    EventInfo[1] = 0;
    EventInfo[2] = Events;
    EventInfo.Write("EventInfo");
    outFile->Close();
    
    //zeroing bins not in data
    tempPions->SetBinContent(17,0);
    tempKaons->SetBinContent(9,0);
    tempKaons->SetBinContent(13,0);
    tempProtons->SetBinContent(7,0);
    tempProtons->SetBinContent(10,0);

    //matching overflow handling in jet analysis
    double newfinaljetbin = HistTempJet->GetBinContent(NpTJetBin)+HistTempJet->GetBinContent(NpTJetBin+1);
	HistTempJet->SetBinContent(NpTJetBin,newfinaljetbin);
    double newfinaldijetbin = HistTempdiJet->GetBinContent(NpTdiJetBin)+HistTempdiJet->GetBinContent(NpTdiJetBin+1);
	HistTempdiJet->SetBinContent(NpTdiJetBin,newfinaldijetbin);

    //Histogram for ration of reco hadrons
    HistRecoHadron->Divide(HistTempSingleHadron);

    //ee Jets debugging
    //cout << integral(HistTempJet) << endl;

    //create root file for total plots
    TFile* totalroot = new TFile( "totals.root", "RECREATE");
    HistTempJet->Write("jets");
    HistTempdiJet->Write("dijets");
    HistTempSingleHadron->Write("hadrons");
    tempPions->Write("pions");
    tempKaons->Write("kaons");
    tempProtons->Write("protons");
    HistThrust->Write("thrust");
    HistMultiplicity->Write("multiplicity");
    HistRecoHadron->Write("Reco Hadron Ratio");

    //ratio plots for comparison
    TFile alephfile("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/data/aleph.root");
    TDirectory* thrustdir = (TDirectory*)alephfile.Get("Table 3");
    TDirectory* multdir = (TDirectory*)alephfile.Get("Table 18");
    TGraphErrors* thrustgraph = (TGraphErrors*) thrustdir->Get("Graph1D_y1");
    TGraphErrors* multgraph = (TGraphErrors*) multdir->Get("Graph1D_y1");
    ratioPlot(thrustgraph, HistThrust, "Thrust", "thrust", "P", false, false);
    ratioPlot(multgraph, HistMultiplicity, "Multiplicity", "N Charged", "P(N)", false, false);

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    
    //debugging
    //for(int i = 0; i < NpTJetBin; i++) cout << DifferentialJetTotal[i]*1000000 << " " << DifferentialJetTotalErrors[i]*1000000 << endl;
    //for(int i = 0; i < NpTSingleHadronBin; i++) cout << DifferentialHadronTotal[i] << " " << DifferentialHadronTotalErrors[i]*1000000 << endl;
    //for(int i = 0; i < NpTHardBin; i++) cout << pTHardBinSingleHadronBinError << endl;

    return 0;
}