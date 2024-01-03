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
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    
    //Reading ptHat bins from list of made directories
    vector<vector<string>> tempvec = getDatBounds("./dat");
    vector<string> pTHatMin = tempvec[0];
    vector<string> pTHatMax = tempvec[1];
    if(stoi(getEcm()) != 2760){
        cout << "Center of mass energy does not match." << endl;
        return 0;
    }
    //for(int i = 0; i < pTHatMin.size(); i++) cout << pTHatMin[i] << endl; //debugging line
    
    //char pTHatMin[][10] = {"0","10","12.5","15","17.5","20","22.5","25","27.5","30","32.5","35","37.5","40","42.5","45","50","55","60","65","70","75","80","90","100","110","120","130","140","150","160","170","180","190","200","210","230","250","270","290","310","330","350","400","450","500","550","600"};
    //char pTHatMax[][10] = {"10","12.5","15","17.5","20","22.5","25","27.5","30","32.5","35","37.5","40","42.5","45","50","55","60","65","70","75","80","90","100","110","120","130","140","150","160","170","180","190","200","210","230","250","270","290","310","330","350","400","450","500","550","600","1000"};
    int NpTHardBin = pTHatMin.size();
    double bin14error[NpTHardBin];
    
    //By default only look at one bin
    //int pTHatMin[1]={110};
    //int pTHatMax[1]={120};
    double DetectorEtaCut= 2.6;
    
    //Variables for jet spectrum
    //(For CMS at 2760 GeV with jet radius R=2,0.3, or 0.4, |eta_jet|<2.0)
    double JetpTBin[] = {64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300}; //in GeV
    int NpTJetBin = sizeof(JetpTBin)/sizeof(JetpTBin[0])-1;     double JetpTMin = 10; //in GeV
    
    vector<int> eventCount;
    double JetEtaCut =2.0;
    double JetRadius = 0.3;
    double dNdpTCountJet[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]
    double pTHardBinJetBinError[NpTHardBin][NpTJetBin];
    double TotalDifferentialJetCrossSection[NpTJetBin] = {0};
    double TotalDifferentialJetCrossSectionError[NpTJetBin] = {0};
    
    //Variables for single hadron spectrum
    double SingleHadronpTBin[] = {0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.5, 1.8, 2.1, 2.4, 3.6, 4.8, 6.0, 7.2, 10.8, 14.4, 21.6, 28.8, 38.4, 48.0, 67.2, 86.4, 112.2};
    int NpTSingleHadronBin = sizeof(SingleHadronpTBin)/sizeof(SingleHadronpTBin[0])-1;

    double pionpTBin[] = {0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
    int NpTpionBin = sizeof(pionpTBin)/sizeof(pionpTBin[0])-1;

    double kaonpTBin[] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
    int NpTkaonBin = sizeof(kaonpTBin)/sizeof(kaonpTBin[0])-1;

    double protonpTBin[] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
    int NpTprotonBin = sizeof(protonpTBin)/sizeof(protonpTBin[0])-1;

    double SingleHadronEtaCut = 1.0;
    double idHadronEtaCut = 0.8;
    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    long double TotalDifferentialSingleHadronYield[NpTSingleHadronBin] = {0};
    long double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin] = {0};
   
    // for jet substructure
    double JetpTCut = 100.0;
    //Variables for jet shape
    double JetShaperBin[7] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
    int NrJetShapeBin = 6;
    //Variables for  jet fragmentation function
    double JetFFzBin[11] = {0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1};
    int NzJetFFBin = 10;
    //Variables for jet mass
    double JetMassBin[13] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
    int NJetMassBin = 12;
    
    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistTempJetShape = new TH1D("JetShapeBin", "Jet Shape", NrJetShapeBin, JetShaperBin); //CountVsr for Jet Shape
    TH1D *HistTempJetFF = new TH1D("JetFragmentaionFunctionBin", "Jet Fragmentaion Function", NzJetFFBin, JetFFzBin);
    TH1D *HistTempJetMass = new TH1D("JetMassBin", "Jet Mass", NJetMassBin, JetMassBin);
	
    TH1D *HistTotalHadron = new TH1D("HadronSpectrumBin", "Combined Hadron pT Spectrum", NpTSingleHadronBin, SingleHadronpTBin); //Total hist for hadrons
	TH1D *HistTotalJet = new TH1D("JetSpectrumBin", "Combined Jet pT Spectrum", NpTJetBin, JetpTBin); //Total hist for jets
    TH1D *HistTotalJetShape = new TH1D("JetShape", "Jet Shape", NrJetShapeBin, JetShaperBin); //Total hist for jet shape
    TH1D *HistTotalJet2 = new TH1D("JetSpectrumBin", "Combined Jet pT Spectrum 0.2 R", NpTJetBin, JetpTBin); //Total hist for jets
	TH1D *HistTotalJet3 = new TH1D("JetSpectrumBin", "Combined Jet pT Spectrum 0.4 R", NpTJetBin, JetpTBin); //Total hist for jets
    TH1D *HistTotalPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin); //identified hadrons hists
    TH1D *HistTotalKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", NpTkaonBin, kaonpTBin);
    TH1D *HistTotalProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protonpTBin);
    HistTotalJet->SetName("Combined Jet pT Spectrum");
    HistTotalHadron->SetName("Combined Hadron pT Spectrum");

    //Running totals for total spectra
    double DifferentialJetTotal[NpTJetBin] = {0};
    double DifferentialJetTotalErrors[NpTJetBin] = {0};
    double DifferentialHadronTotal[NpTSingleHadronBin] = {0};
    double DifferentialHadronTotalErrors[NpTSingleHadronBin] = {0};
    double JetShapeTotals[NrJetShapeBin] = {0,0,0,0,0,0};
    
    std::vector <fjcore::PseudoJet> fjInputs;
    std::vector <int> chargeList;
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);
    fjcore::JetDefinition jetDef2(fjcore::antikt_algorithm, 0.2);
    fjcore::JetDefinition jetDef3(fjcore::antikt_algorithm, 0.4);
    
    Pythia8::Pythia pythia;//("",false);
    
    //specific bins to examine
    NpTHardBin = 2;
    pTHatMin[0] = "0"; pTHatMin[1] = "0";
    pTHatMax[0] = "5"; pTHatMax[1] = "5";
    string HadronFiles[] = {"dat/PP_Bin0_5new.dat","dat/PP_Bin0_5new.dat"};
	
	//get list of cross sections
    vector<vector<double>> xsecout = getXsecs(pTHatMin,pTHatMax);
	vector<double> xsecList = xsecout[0];
    vector<double> xsecErrorList = xsecout[0];
    double xsectotal = 0;
    for(int k = 0; k<NpTHardBin; k++) xsectotal += xsecList[k];
    //for(int k = 0; k<NpTHardBin; k++) cout << pTHatMin[k] << " " << pTHatMax[k] << " " << xsecList[k]*100000 << endl; //debugging line

    //graph declaration for adding hadron spectra
    TMultiGraph* hadronComp = new TMultiGraph();
    TGraph* hadronComponents[NpTHardBin];
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char pTBinString[100];
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFiles[k].c_str());
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass;
        int Events =0;
        int TriggeredJetNumber=0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s_%i.root",pTHatMin[k].c_str(),pTHatMax[k].c_str(),k);
        TFile* outFile = new TFile( outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];
        
        HistTempJet->Reset();
        sprintf(HistName,"CountVspTJetSpectrumBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempJet->SetName(HistName);
        
        HistTempSingleHadron->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempSingleHadron->SetName(HistName);
        
        HistTempJetShape->Reset();
        sprintf(HistName,"CountVsrJetShapeBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempJetShape->SetName(HistName);
        HistTempJetShape->Sumw2();
        
        HistTempJetFF->Reset();
        sprintf(HistName,"CountVszJetFFBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempJetFF->SetName(HistName);
        HistTempJetFF->Sumw2();
        
        HistTempJetMass->Reset();
        sprintf(HistName,"CountVsJetMassBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempJetMass->SetName(HistName);
        HistTempJetMass->Sumw2();
        
        fjInputs.resize(0);
        chargeList.resize(0);
        
        //temp hists for jets
        TH1D *HistTempJet2 = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin);
        TH1D *HistTempJet3 = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin);

        //temp hists for identified hadrons
        TH1D *tempPions = new TH1D("Pion Spectrum", "Pion Spectrum pT", NpTpionBin, pionpTBin);
        TH1D *tempKaons = new TH1D("Kaon Spectrum", "Kaon Spectrum pT", NpTkaonBin, kaonpTBin);
        TH1D *tempProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protonpTBin);

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        double hadskept[] = {0,0,0,0}; double hadslost[] ={0,0,0,0}; double scale [4]; //handling scaling for lost hadrons in blending
        vector <fjcore::PseudoJet> UnSortedJets, SortedJets, UnSortedJets2, SortedJets2, UnSortedJets3, SortedJets3, constituents;
        
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
                double PT = TMath::Sqrt( (Px*Px) + (Py*Py));
                
                // Add this particle for Jet spectrum
                if( SN==0 && fjInputs.size()>0 ){
                    
                    // List first few FastJet jets and some info about them.
                    /* if (nListJets)
                    {
                        cout << "\n --------  FastJet jets, R = " << JetRadius << "anti-kt for pTHatBin = "<<k
                        << "  --------------------------------------------------\n\n "
                        << "  i         pT        y     eta      phi  " << endl;
                        for (int i = 0; i < int(SortedJets.size()); ++i)
                        {
                            vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();
                            cout << setw(4) << i << fixed << setprecision(3) << setw(11)
                            << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap()
                            << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << endl;
                        }
                        cout << "\n --------  End FastJet Listing  ------------------"
                        << "---------------------------------" << endl;
                    } */
                    
                    //Alternate Radius Calculations, 0.2 first
                    fjcore::ClusterSequence clustSeq2(fjInputs, jetDef2);
                    UnSortedJets2 = clustSeq2.inclusive_jets(JetpTMin);
                    SortedJets2 = sorted_by_pt(UnSortedJets2);
                    int pFast2 = SortedJets2.size();
                    for (int i = 0; i < pFast2; ++i){
                        if(-JetEtaCut < SortedJets2[i].eta() && SortedJets2[i].eta()< JetEtaCut){
                            if(SortedJets2[i].perp() < stod(pTHatMax[k])*1.1) HistTempJet2->Fill(SortedJets2[i].perp());
                            else HistTempJet2->Fill(stod(pTHatMax[k])); //filtering out anomalous high energy events
                            //cout << "jet rad 0.2 " << SortedJets2[i].perp() << endl;
                        }
                    }

                    //jet radius 0.4
                    fjcore::ClusterSequence clustSeq3(fjInputs, jetDef3);
                    UnSortedJets3 = clustSeq3.inclusive_jets(JetpTMin);
                    SortedJets3 = sorted_by_pt(UnSortedJets3);
                    int pFast3 = SortedJets3.size();
                    for (int i = 0; i < pFast3; ++i){
                        if(-JetEtaCut < SortedJets3[i].eta() && SortedJets3[i].eta()< JetEtaCut){
                            if(SortedJets3[i].perp() < stod(pTHatMax[k])*1.1) HistTempJet3->Fill(SortedJets3[i].perp());
                            else HistTempJet3->Fill(stod(pTHatMax[k])); //filtering out anomalous high energy events
                            //cout << "jet rad 0.4" << endl;
                        }
                    }

                    
                    // Run Fastjet algorithm and sort jets in pT order.
                    fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
                    UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
                    SortedJets    = sorted_by_pt(UnSortedJets);

                    int pFast = SortedJets.size();
                    for (int i = 0; i < pFast; ++i)
                    {
                        if(-JetEtaCut < SortedJets[i].eta() && SortedJets[i].eta()< JetEtaCut )
                        {
                            if(SortedJets[i].perp() < stod(pTHatMax[k])*1.1) HistTempJet->Fill( SortedJets[i].perp() );
                            else HistTempJet->Fill(stod(pTHatMax[k])); //filtering out anomalous high energy events

                            if(SortedJets[i].perp() >= JetpTCut && SortedJets[i].perp() < stod(pTHatMax[k])*1.1){

                                TriggeredJetNumber ++;
                                constituents = SortedJets[i].constituents();

                                //Jet Shape---
                                for( int j = 0; j < constituents.size(); j++ ){
                                    
                                    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    HistTempJetShape->Fill( delta_r, constituents[j].perp() );
                                }
                                
                                //Jet Fragmentation Function---
                                for( int j = 0; j < fjInputs.size(); j++ ){
                                    double delta_eta = fjInputs[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(fjInputs[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    if( fabs(chargeList[j]) > 0.01 && delta_r <= JetRadius){//charged particle in jet cone
                                        double z_jet = fjInputs[j].perp()/SortedJets[i].perp();
                                        HistTempJetFF->Fill( z_jet );
                                    }
                                }

                                //Jet Mass---
                                double jet_e = 0.0, jet_px = 0.0, jet_py = 0.0, jet_pz = 0.0;
                                for( int j = 0; j < constituents.size(); j++ ){
                                    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    if( delta_r <= JetRadius){// "ALL" particle in jet cone
                                        jet_e += constituents[j].e();
                                        jet_px += constituents[j].px();
                                        jet_py += constituents[j].py();
                                        jet_pz += constituents[j].pz();
                                    }
                                }
                                double jet_mass
                                = TMath::Sqrt(jet_e*jet_e - jet_px*jet_px - jet_py*jet_py - jet_pz*jet_pz);
                                HistTempJetMass->Fill( jet_mass );

                            } //commented out to reduce run time per event
                        }
                    }
                    fjInputs.resize(0);
                    //cout<<"Found a Jet \t "<<pTBinString<<"\t NetJetevents is \t"<<NetJetEvents<<endl;
                    if(  fabs(Eta) < DetectorEtaCut &&  PT>0.01  && PID!=12 && PID!=14 && PID!=16 && PID!=18){
                        fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                        chargeList.push_back( pythia.particleData.charge( PID ) );
                    }
                }else{
                    //cout<<EventLabel << " "<<PID <<" 1"<< E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;
                    if( fabs(Eta) < DetectorEtaCut && PT>0.01  &&  PID!=12 && PID!=14 && PID!=16 && PID!=18 ){
                        fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                        chargeList.push_back( pythia.particleData.charge( PID ) );
                    }
                }
                
                
                // Add this particle into SingleHadron spectrum
                if(fabs(Eta) < SingleHadronEtaCut && PT>0.01  && fabs(PID)>100 &&  pythia.particleData.charge(PID)!=0){
                    //cout<<PT<<" PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
                    
                    //blending weights
                    double strength = 1;//blending(k,PT,stod(pTHatMax[k]),1,stod(pTHatMax[0]));
                    hadskept[0] += strength; hadslost[0] += 1-strength;

                    HistTempSingleHadron->Fill(PT, strength); 
                    strength *= sqrt(pow(mass,2) + pow(PT*cosh(Eta),2))/(PT*cosh(Eta)); //corrective factor for rapidity graph, must be done after total hadrons graph
                    if(fabs(Eta) < idHadronEtaCut){
                        if(abs(PID) == 211) {tempPions->Fill(PT, strength); hadskept[1] += strength; hadslost[1] += 1-strength;}
                        if(abs(PID) == 321) {tempKaons->Fill(PT, strength); hadskept[2] += strength; hadslost[2] += 1-strength;}
                        if(abs(PID) == 2212) {tempProtons->Fill(PT, strength); hadskept[3] += strength; hadslost[3] += 1-strength;}
                    } 
                }
            }
        }

        //xsec and event count handling
        eventCount.push_back(Events);
        double HardCrossSection = xsecList[k];
        double HardCrossSectionError = xsecErrorList[k];
        for(int i = 0; i < 4; i++){
            scale[i] = 1 + (hadslost[i]/hadskept[i]);
            hadslost[i] = 0; hadskept[i] = 0;
        }
        
        // end of reading cross section
        
        //For jet spectrum, weighted by cross section and combined through multiple pTHatBins
        for(int j=0;j<NpTJetBin;j++)
        {
            dNdpTCountJet[k][j]= HistTempJet->GetBinContent(j+1);
            if(dNdpTCountJet[k][j] > 0.0)
            {
                pTHardBinJetBinError[k][j] = (dNdpTCountJet[k][j]*HardCrossSection/(Events*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut))*TMath::Sqrt( (1/dNdpTCountJet[k][j]) + TMath::Power(HardCrossSectionError/HardCrossSection,2.0));
                cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountJet[k][j]*HardCrossSection)/(Events*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut)<<endl;
            }
            else
            {
                pTHardBinJetBinError[k][j] = 0.0;
                cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
            }
        }

        //filtering out low bin events in the hadron spectra, smoothing
        for(int i = 0; i < NpTSingleHadronBin; i++)
            if(HistTempSingleHadron->GetBinContent(i+1) < 5) HistTempSingleHadron->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTpionBin; i++)
            //if(tempPions->GetBinContent(i+1) < 2) tempPions->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTkaonBin; i++)
            //if(tempKaons->GetBinContent(i+1) < 2) tempKaons->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTpionBin; i++)
            //if(tempProtons->GetBinContent(i+1) < 2) tempProtons->SetBinContent(i+1,0);
        
        
        //For single Hadron spectrum
        for(int j=0;j<NpTSingleHadronBin;j++)
        {
            dNdpTCountSingleHadron[k][j]= HistTempSingleHadron->GetBinContent(j+1);
            if(dNdpTCountSingleHadron[k][j]> 0.0)
            {
                pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection/(Events*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError/HardCrossSection,2.0));
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection)/(Events*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
            }
            else
            {
                pTHardBinSingleHadronBinError[k][j] = 0.0;
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
            }
        }
        
        //Write histogram into a root file
        HistTempJet->Write();
        HistTempJetShape->Write();
        HistTempSingleHadron->Write();
        HistTempJetFF->Write();
        HistTempJetMass->Write();
        tempPions->Write();
        tempKaons->Write();
        tempProtons->Write();
        
        //add to totals histograms
        HistTotalHadron->Add(HistTempSingleHadron,scale[0]*HardCrossSection/Events);
        HistTotalPions->Add(tempPions,scale[1]*HardCrossSection/Events);
        HistTotalKaons->Add(tempKaons,scale[2]*HardCrossSection/Events);
        HistTotalProtons->Add(tempProtons,scale[3]*HardCrossSection/Events);
		HistTotalJet->Add(HistTempJet,HardCrossSection); 
        HistTotalJet2->Add(HistTempJet2,HardCrossSection/Events);
        HistTotalJet3->Add(HistTempJet3,HardCrossSection/Events);
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");
        
        TVector TriggeredJetInfo(2);
        TriggeredJetInfo[0] = JetpTCut;
        TriggeredJetInfo[1] = TriggeredJetNumber;
        TriggeredJetInfo.Write("TriggeredJetInfo");
        
        //Plots for jet spectrum
        double DifferentialJetCrossSection[NpTJetBin],DifferentialJetCrossSectionError[NpTJetBin],JetpT[NpTJetBin],JetpTError[NpTJetBin];
        TGraphErrors * GEJet;
        
        cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
        for(int j=0; j<NpTJetBin;j++){
            DifferentialJetCrossSection[j] = (dNdpTCountJet[k][j]*HardCrossSection)/(Events*((JetpTBin[j+1]-JetpTBin[j]))*2.0*JetEtaCut);
            DifferentialJetCrossSectionError[j] = pTHardBinJetBinError[k][j];
            JetpT[j] = (JetpTBin[j] + JetpTBin[j+1])/2.0;
            JetpTError[j] = (JetpTBin[j+1] - JetpTBin[j])/2.0;
            cout<<JetpT[j]<<"\t"<<DifferentialJetCrossSection[j]<<"\t"<<DifferentialJetCrossSectionError[j]<<"\t"<<dNdpTCountJet[k][j]<<"\t"<<HardCrossSection<<endl;
            DifferentialJetTotal[j] += DifferentialJetCrossSection[j];
            DifferentialJetTotalErrors[j] += DifferentialJetCrossSectionError[j]*DifferentialJetCrossSectionError[j];
        }
        
        GEJet = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetCrossSection,JetpTError,DifferentialJetCrossSectionError);
        char MyGraphName[100];
        sprintf(MyGraphName,"DifferentialJetCrossSectionBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        GEJet->SetNameTitle(MyGraphName);
        GEJet->Write();
        
        // For charged Hadron spectrum
        double DifferentialSingleHadronYield[NpTSingleHadronBin],DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
        TGraphErrors * GESingleHadron;
        
        cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
        for(int j=0; j<NpTSingleHadronBin;j++){
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection)/(Events*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*xsectotal*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*2.0*SingleHadronEtaCut); //hadron scaling, changed to scale with total cross section from pT hard bins
            DifferentialSingleHadronYieldError[j] = pTHardBinSingleHadronBinError[k][j];
            SingleHadronpT[j] = (SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1]-SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYield[j]<<"\t"<<DifferentialSingleHadronYieldError[j]<<endl;
            DifferentialHadronTotal[j] += DifferentialSingleHadronYield[j];
            DifferentialHadronTotalErrors[j] += DifferentialSingleHadronYieldError[j]*DifferentialSingleHadronYieldError[j];
        }
        
        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        GESingleHadron->SetNameTitle(MyGraphName2);
        GESingleHadron->Write();
		
		//Save jet pT hist as a png for convenience
		/*char OutJetHistName[100];
		sprintf(OutJetHistName,"plots/JetSpectrumBin%s_%s.png",pTHatMin[k].c_str(),pTHatMax[k].c_str());
		TCanvas *cJet = new TCanvas();
        cJet->SetLogy();
		//GEJet->Draw();
		//cJet->Print(OutJetHistName);
        if (cJet) {cJet->Close();}*/

        //Save hadron pT hist as a png for convenience
        hadronComponents[k] = (TGraph*)GESingleHadron->Clone();
        hadronComponents[k]->SetLineColor(k+2);
		/*char OutHadronHistName[100];
		sprintf(OutHadronHistName,"plots/HadronSpectrumBin%s_%s.png",pTHatMin[k].c_str(),pTHatMax[k].c_str());
		TCanvas *cHadron = new TCanvas();
        cHadron->SetLogy();
        cHadron->SetLogx();
		GESingleHadron->Draw();
		cHadron->Print(OutHadronHistName);
        if (cHadron) {cHadron->Close();}*/
        
        // for jet shape
        HistTempJetShape->Scale( (1.0/TriggeredJetNumber), "width" );
        double ErrorNormJetShape;
        double NormJetShape
        = HistTempJetShape->IntegralAndError( 1, NrJetShapeBin,
                                              ErrorNormJetShape,
                                              "width" );
        TH1D *Norm = (TH1D*)HistTempJetShape->Clone("Normalization");
        Norm->Sumw2();Norm->SetTitle("Normalization");
        int nbins = Norm->GetSize();
        for( int i=1; i < nbins-1; i++){
            Norm->SetBinContent(i, NormJetShape);
            Norm->SetBinError(i, ErrorNormJetShape);
        }
        HistTempJetShape->Divide(Norm);
        delete Norm;

        double JetShape[NrJetShapeBin],JetShapeError[NrJetShapeBin],JetShapeR[NrJetShapeBin],JetShapeRError[NrJetShapeBin];
        TGraphErrors * GEJetShape;
        
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Shape is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NrJetShapeBin;j++){
            JetShape[j] = HistTempJetShape->GetBinContent(j+1);
            double jetshapetemp = JetShape[j] * HardCrossSection;
            if(isnan(jetshapetemp) == 0) JetShapeTotals[j] += jetshapetemp; 
            JetShapeError[j] = HistTempJetShape->GetBinError(j+1);
            JetShapeR[j] = (JetShaperBin[j]+JetShaperBin[j+1])/2.0;
            JetShapeRError[j] = (JetShaperBin[j+1]-JetShaperBin[j])/2.0;
            cout<<JetShapeR[j]<<"\t"<<JetShape[j]<<"\t"<<JetShapeError[j]<<endl;
        }
        GEJetShape = new TGraphErrors(NrJetShapeBin,JetShapeR,JetShape,JetShapeRError,JetShapeError);
        char MyGraphName3[100];
        sprintf(MyGraphName3,"JetShapeBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        GEJetShape->SetName(MyGraphName3);
        GEJetShape->Write();
        
        // for jet fragmentation function
        HistTempJetFF->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetFF[NzJetFFBin], JetFFError[NzJetFFBin], JetZ[NzJetFFBin], JetZError[NzJetFFBin];
        TGraphErrors * GEJetFF;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Fragmentation Function is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NzJetFFBin;j++){
            JetFF[j] = HistTempJetFF->GetBinContent(j+1);
            JetFFError[j] = HistTempJetFF->GetBinError(j+1);
            JetZ[j] = (JetFFzBin[j]+JetFFzBin[j+1])/2.0;
            JetZError[j] = (JetFFzBin[j+1]-JetFFzBin[j])/2.0;
            cout<<JetZ[j]<<"\t"<<JetFF[j]<<"\t"<<JetFFError[j]<<endl;
        }
        GEJetFF = new TGraphErrors(NzJetFFBin,JetZ,JetFF,JetZError,JetFFError);
        char MyGraphName4[100];
        sprintf(MyGraphName4,"JetFFBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        GEJetFF->SetName(MyGraphName4);
        GEJetFF->Write();
        
        // for jet mass
        HistTempJetMass->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetMassDist[NJetMassBin], JetMassDistError[NJetMassBin], JetMass[NJetMassBin], JetMassError[NJetMassBin];
        TGraphErrors * GEJetMass;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Mass Distribution is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NJetMassBin;j++){
            JetMassDist[j] = HistTempJetMass->GetBinContent(j+1);
            JetMassDistError[j] = HistTempJetMass->GetBinError(j+1);
            JetMass[j] = (JetMassBin[j]+JetMassBin[j+1])/2.0;
            JetMassError[j] = (JetMassBin[j+1]-JetMassBin[j])/2.0;
            cout<<JetMass[j]<<"\t"<<JetMassDist[j]<<"\t"<<JetMassDistError[j]<<endl;
        }
        GEJetMass = new TGraphErrors(NJetMassBin,JetMass,JetMassDist,JetMassError,JetMassDistError);
        char MyGraphName5[100];
        sprintf(MyGraphName5,"JetMassDistBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        GEJetMass->SetName(MyGraphName5);
        GEJetMass->Write();

        delete GEJet; delete GESingleHadron; delete GEJetShape; delete GEJetFF; delete GEJetMass;
        delete outFile;
    } //k-loop ends here (pTHatBin loop)
    
    delete HistTempJet; delete HistTempSingleHadron;
    delete HistTempJetShape; delete HistTempJetFF; delete HistTempJetMass;

    //create histogram for ratio plot
    TH1D *HistFinalJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistFinalSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron

    //setting bin centers and widths, rescaling appropriately
    double JetpT[NpTJetBin]; double JetpTError[NpTJetBin];
    for(int i = 0; i < NpTJetBin; i++) {
        JetpT[i] = (JetpTBin[i] + JetpTBin[i+1])/2.0;
        JetpTError[i] = (JetpTBin[i+1] - JetpTBin[i])/2.0;
        DifferentialJetTotalErrors[i] = sqrt(DifferentialJetTotalErrors[i]);
        if(DifferentialJetTotalErrors[i] > DifferentialJetTotal[i]) DifferentialJetTotalErrors[i] = DifferentialJetTotal[i]*0.99; //correcting for negatives
        HistFinalJet->SetBinContent(i+1,DifferentialJetTotal[i]);//for ratio plot
        HistFinalJet->SetBinError(i+1,DifferentialJetTotalErrors[i]);
    }

    double SingleHadronpT[NpTSingleHadronBin]; double SingleHadronpTError[NpTHardBin];
    for(int i = 0; i < NpTSingleHadronBin; i++){
        SingleHadronpT[i] = (SingleHadronpTBin[i]+SingleHadronpTBin[i+1])/2.0;
        SingleHadronpTError[i] = (SingleHadronpTBin[i+1]-SingleHadronpTBin[i])/2.0;
        DifferentialHadronTotalErrors[i] = sqrt(DifferentialHadronTotalErrors[i]);
        if(DifferentialHadronTotalErrors[i] > DifferentialHadronTotal[i])  DifferentialHadronTotalErrors[i] = DifferentialHadronTotal[i]*0.99; //correcting for negatives
        HistFinalSingleHadron->SetBinContent(i+1, DifferentialHadronTotal[i]); //for ratio plot
        HistFinalSingleHadron->SetBinError(i+1, DifferentialHadronTotalErrors[i]);
    }

    //Scaling totals by global factors and the identified pions by bin centers
    HistTotalJet2->Scale(1/(2.0*JetEtaCut),"width");
    HistTotalJet3->Scale(1/(2.0*JetEtaCut),"width");
    //scaleBins(HistTotalHadron,1/(xsectotal*2*M_PI*2.0*SingleHadronEtaCut));
    //scaleBins(HistTotalPions,1/(xsectotal*2*M_PI*2.0*idHadronEtaCut));
    //scaleBins(HistTotalKaons,1/(xsectotal*2*M_PI*2.0*idHadronEtaCut));
    //scaleBins(HistTotalProtons,1/(xsectotal*2*M_PI*2.0*idHadronEtaCut));

    //jet shape
    double JetShapeNorm = 0;
    for(int i = 0; i < NrJetShapeBin; i++){
        JetShapeNorm += JetShapeTotals[i]*(JetShaperBin[1]-JetShaperBin[0]);
        HistTotalJetShape->SetBinContent(i+1, JetShapeTotals[i]);
    }
    HistTotalJetShape->Scale(1/JetShapeNorm);

    //jets
    //jet data graph
    TMultiGraph *GEJetTotal = new TMultiGraph();
    TFile jet_file("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/data/JetData.root");
    TDirectory* jetdir = (TDirectory*)jet_file.Get("Table 4");
    TGraphErrors* jetData = (TGraphErrors*) jetdir->Get("Graph1D_y2");
    TH1D* jetDataHistOriginal = (TH1D*) jetdir->Get("Hist1D_y2");
    TH1D* jetDataHist = (TH1D*) jetDataHistOriginal->Clone(); //cloning to avoid messing up original
    jetDataHist->Scale(1e-6); //data listed in mb
    jetData->SetMarkerStyle(kCircle);
    jetData->SetMarkerColor(kRed);
    jetData->SetLineColor(kRed);

    double jetDataHistErrors[] = {0.35,0.15,0.066,0.065,0.026,0.0099,0.0038,0.0016,0.00065,0.00028,0.00013,5.7e-5};
    for(int i = 0; i < jetDataHist->GetNbinsX(); i++) jetDataHist->SetBinError(i+1,jetDataHistErrors[i]*1e-6);
    //ratioPlot(jetDataHist,HistFinalJet,"Jet r3 Differential Cross sections (mb/GeV/c)");

    //r2 graph
    TH1D* jetDataHistOriginal2 = (TH1D*) jetdir->Get("Hist1D_y1");
    TH1D* jetDataHist2 = (TH1D*) jetDataHistOriginal2->Clone(); //cloning to avoid messing up original
    jetDataHist2->Scale(1e-6); //data listed in mb

    double jetDataHistErrors2[] = {0.24,0.1,0.048,0.048,0.02,0.0074,0.0028,0.0012,0.00054,0.00024,0.00011,4.6e-5};
    for(int i = 0; i < jetDataHist2->GetNbinsX(); i++) jetDataHist2->SetBinError(i+1,jetDataHistErrors2[i]*1e-6);
    //ratioPlot(jetDataHist2,HistTotalJet2,"Jet r2 Differential Cross sections (mb/GeV/c)");

    //r4 graph
    TH1D* jetDataHistOriginal4 = (TH1D*) jetdir->Get("Hist1D_y3");
    TH1D* jetDataHist4 = (TH1D*) jetDataHistOriginal4->Clone(); //cloning to avoid messing up original
    jetDataHist4->Scale(1e-6); //data listed in mb

    double jetDataHistErrors4[] = {0.45,0.18,0.083,0.082,0.033,0.012,0.0045,0.0019,0.00078,0.00034,0.00016,7e-5};
    for(int i = 0; i < jetDataHist4->GetNbinsX(); i++) jetDataHist4->SetBinError(i+1,jetDataHistErrors4[i]*1e-6);
    //ratioPlot(jetDataHist4,HistTotalJet3,"Jet r4 Differential Cross sections (mb/GeV/c)");

    //jet jetscape graph
    TGraphErrors *GEJetPrediction;
    GEJetPrediction = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetTotal,JetpTError,DifferentialJetTotalErrors);
    GEJetPrediction->SetMarkerStyle(kCircle);
    GEJetPrediction->SetLineColor(kBlue);
    GEJetPrediction->SetMarkerColor(kBlue);

    //jet total
    GEJetTotal->Add(GEJetPrediction);
    GEJetTotal->Add(jetData);
    GEJetTotal->SetName("Total Jet Spectra");
    GEJetTotal->SetTitle("Total Jet Spectra");
    GEJetTotal->GetXaxis()->SetTitle("Jet pT (GeV)");
    GEJetTotal->GetYaxis()->SetTitle("Differential Cross-Section (c*mb/GeV)");

    //jet canvas and saving
    TCanvas *cTotalJet = new TCanvas();
    cTotalJet->SetLogy();
    TLegend jetleg(.7,.7,.9,.9,"Jet Spectra");
    jetleg.AddEntry(GEJetPrediction,"Jetscape","lep");
    jetleg.AddEntry(jetData,"CMS Data","lep");
	//GEJetTotal->Draw("apl");
    //jetleg.DrawClone("Same");
    auto rpJet = new TRatioPlot(HistFinalJet, jetDataHist);
    //rpJet->Draw("C");
	//cTotalJet->Print("plots/JetTotalSpectrum.png");
    jet_file.Close();

    //hadrons
    //hadron data graph
    TMultiGraph *GEHadronTotal = new TMultiGraph();
    TFile hadron_file("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/data/HadronData.root");
    TDirectory* hadrondir = (TDirectory*)hadron_file.Get("Table 1");
    TGraphErrors* hadronData = (TGraphErrors*) hadrondir->Get("Graph1D_y1");
    TH1D* hadronDataHist = (TH1D*)hadrondir->Get("Hist1D_y1");
    hadronData->SetMarkerStyle(kCircle);
    hadronData->SetMarkerColor(kRed);
    hadronData->SetLineColor(kRed);

    double hadDataHistErrors[] = {0.0589,0.0306,0.0166,0.00958,0.00576,0.00285,0.00122,0.000565,0.00028,6e-5,8.17e-6,1.73e-6,4.69e-7,6.34e-8,7.04e-9,6.85e-10,8.27e-11,1.43e-11,3.26e-12,5.54e-13,9.02e-14,1.59e-14};
    for(int i = 0; i < hadronDataHist->GetNbinsX(); i++) hadronDataHist->SetBinError(i+1,hadDataHistErrors[i]);
    //ratioPlot(hadronDataHist, HistFinalSingleHadron, "Hadron Yields (Gev^-2*c^3)", true);

    //hadron jetscape
    TGraphErrors *GEHadronPrediction;
    GEHadronPrediction = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialHadronTotal,SingleHadronpTError,DifferentialHadronTotalErrors);
    GEHadronPrediction->SetMarkerStyle(kDot);
    GEHadronPrediction->SetMarkerColor(kBlue);
    GEHadronPrediction->SetLineColor(kBlue);
    GEHadronPrediction->GetXaxis()->SetTitle("Hadron pT (GeV)");
    GEHadronPrediction->GetYaxis()->SetTitle("Hadron pT (GeV)");
    for(int i = 0; i < NpTHardBin; i++) hadronComp->Add(hadronComponents[i]);
    hadronComp->Add(GEHadronPrediction);
    hadronComp->GetYaxis()->SetRangeUser(1e-15,10);

    //hadron total
    GEHadronTotal->Add(GEHadronPrediction);
    GEHadronTotal->Add(hadronData);
    GEHadronTotal->SetName("Total Hadron Spectra");
    GEHadronTotal->SetTitle("Total Hadron Spectra");
    GEHadronTotal->GetXaxis()->SetTitle("Hadron pT (GeV)");
    GEHadronTotal->GetYaxis()->SetTitle("Hadron Yields (Gev^-2*c^3)");

    //hadron draw and saving
    TCanvas *cTotalHadron = new TCanvas();
    cTotalHadron->SetLogy();
    cTotalHadron->SetLogx();
    TLegend hadronleg(.7,.7,.9,.9,"Hadron Yields");
    hadronleg.AddEntry(GEHadronPrediction,"Jetscape","lep");
    hadronleg.AddEntry(hadronData,"CMS Data","lep");
    for(int i = 0; i < NpTHardBin; i++) hadronleg.AddEntry(hadronComponents[i],pTHatMax[i].c_str(),"lep");
	hadronComp->Draw("alX");
    hadronleg.DrawClone("Same");
	cTotalHadron->Print("plots/HadronTotalSpectrumtest.png");
    hadron_file.Close();
	
    //create root file for total plots
    TFile* totalroot = new TFile( "root/test.root", "RECREATE");
    GEJetPrediction->Write("jets graph");
    GEHadronPrediction->Write("hadrons graph");
    HistFinalJet->Write("jet radius 0.3");
    HistTotalHadron->Write("hadrons");
    HistTotalJetShape->Write("jet shape");
    HistTotalJet2->Write("jet radius 0.2");
    HistTotalJet3->Write("jet radius 0.4");
    HistTotalPions->Write("raw pions"); smoothBins(HistTotalPions); HistTotalPions->Smooth(); HistTotalPions->Write("identified pions");
    HistTotalKaons->Write("raw kaons"); smoothBins(HistTotalKaons); HistTotalKaons->Smooth(); HistTotalKaons->Write("identified kaons");
    HistTotalProtons->Write("raw protons"); smoothBins(HistTotalProtons); HistTotalProtons->Smooth(); HistTotalProtons->Write("identified protons");
    hadronComp->Write("hadron components");

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
