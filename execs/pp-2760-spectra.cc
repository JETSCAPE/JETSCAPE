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
#include <TDirectory.h>
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
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");
    
    //Reading ptHat bins from list of made directories
    vector<vector<string>> tempvec = getDatBounds("./dat");
    vector<string> pTHatMin = tempvec[0];
    vector<string> pTHatMax = tempvec[1];
    int NpTHardBin = pTHatMin.size();
    double bin14error[NpTHardBin];
    //for(int i = 0; i < pTHatMin.size(); i++) cout << pTHatMin[i] << endl; //debugging line
    
    
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

    //ID Had variables
    TFile idhadron_file("/data/rjfgroup/rjf01/cameron.parker/data/LHC-ID-hads.root");
    TDirectory* piondir = (TDirectory*)idhadron_file.Get("Table 1");
    TH1D* piondata = (TH1D*) piondir->Get("Hist1D_y3");
    TGraphErrors* piongraph = (TGraphErrors*) piondir->Get("Graph1D_y3");
    int NpTpionBin = piondata->GetNbinsX();
    
    TDirectory* kaondir = (TDirectory*)idhadron_file.Get("Table 2");
    TH1D* kaondata = (TH1D*) kaondir->Get("Hist1D_y3");
    TGraphErrors* kaongraph = (TGraphErrors*) kaondir->Get("Graph1D_y3");
    int NpTkaonBin = kaondata->GetNbinsX();
    
    TDirectory* protondir = (TDirectory*)idhadron_file.Get("Table 3");
    TH1D* protondata = (TH1D*) protondir->Get("Hist1D_y3");
    TGraphErrors* protongraph = (TGraphErrors*) protondir->Get("Graph1D_y3");
    int NpTprotonBin = protondata->GetNbinsX();

    double SingleHadronEtaCut = 1.0;
    double DetectorEtaCut= 2.6;
    double idHadronYCut = 0.8;
    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    long double TotalDifferentialSingleHadronYield[NpTSingleHadronBin] = {0};
    long double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin] = {0};
    double softend = 6.0;
   
    // for jet substructure
    double JetpTCut = 100.0;
    
    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempSingleHadronSoft = new TH1D("SingleHadronSpectrumBinSoft", "Soft Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistTempSingleHadronHard = new TH1D("SingleHadronSpectrumBinHard", "Hard Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
	
    TH1D *HistTotalHadronSoft = new TH1D("HadronSpectrumBinSoft", "Soft Hadron pT Spectrum", NpTSingleHadronBin, SingleHadronpTBin); //Total hist for hadrons
    TH1D *HistTotalHadronHard = new TH1D("HadronSpectrumBinHard", "Hard Hadron pT Spectrum", NpTSingleHadronBin, SingleHadronpTBin); //Total hist for hadrons
	TH1D *HistTotalJet = new TH1D("JetSpectrumBin1", "Combined Jet pT Spectrum", NpTJetBin, JetpTBin); //Total hist for jets
    TH1D *HistTotalJet2 = new TH1D("JetSpectrumBin2", "Combined Jet pT Spectrum 0.2 R", NpTJetBin, JetpTBin); //Total hist for jets
	TH1D *HistTotalJet3 = new TH1D("JetSpectrumBin3", "Combined Jet pT Spectrum 0.4 R", NpTJetBin, JetpTBin); //Total hist for jets
    TH1D *HistTotalPionsSoft = new TH1D("Pion Spectrum Soft", "Pion Spectrum pT", NpTpionBin, piondata->GetXaxis()->GetXbins()->GetArray()); //identified hadrons hists
    TH1D *HistTotalPionsHard = new TH1D("Pion Spectrum Hard", "Pion Spectrum pT", NpTpionBin, piondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalKaonsSoft = new TH1D("Kaon Spectrum Soft", "Kaon Spectrum pT", NpTkaonBin, kaondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalKaonsHard = new TH1D("Kaon Spectrum Hard", "Kaon Spectrum pT", NpTkaonBin, kaondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalProtonsSoft = new TH1D("Proton Spectrum Soft", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalProtonsHard = new TH1D("Proton Spectrum Hard", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistRecoProtons = new TH1D("Proton Spectrum", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());
    HistTotalJet->SetName("Combined Jet pT Spectrum");
    //HistTotalHadron->SetName("Combined Hadron pT Spectrum");

    //Running totals for total spectra
    double DifferentialJetTotal[NpTJetBin] = {0};
    double DifferentialJetTotalErrors[NpTJetBin] = {0};
    double DifferentialHadronTotal[NpTSingleHadronBin] = {0};
    double DifferentialHadronTotalErrors[NpTSingleHadronBin] = {0};
    
    std::vector <fjcore::PseudoJet> fjInputs;
    std::vector <int> chargeList;
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);
    fjcore::JetDefinition jetDef2(fjcore::antikt_algorithm, 0.2);
    fjcore::JetDefinition jetDef3(fjcore::antikt_algorithm, 0.4);
    
    Pythia8::Pythia pythia;//("",false);
	
	//get list of cross sections
    double xsectotal = 62.8; //experimental value: https://arxiv.org/pdf/1208.4968.pdf
    //for(int k = 1; k<NpTHardBin; k++) xsectotal += xsecList[k]; //skip soft bin since it is not added to hard spectra
    //for(int k = 0; k<NpTHardBin; k++) cout << pTHatMin[k] << " " << pTHatMax[k] << " " << xsecList[k]*100000 << endl; //debugging line

    //graph declaration for adding hadron spectra
    TMultiGraph* hadronComp = new TMultiGraph();
    TGraph* hadronComponents[NpTHardBin];
    TDirectory* binfiles[NpTHardBin];
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat.gz", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Y, Phi, pStat, mass;
        int Events =0;
        int TriggeredJetNumber=0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"SpectraBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        binfiles[k] = totalroot->mkdir(outFileName);
        binfiles[k]->cd();
        // Reset for each pTHardBin
        char HistName[100];
        
        HistTempJet->Reset();
        sprintf(HistName,"CountVspTJetSpectrumBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        HistTempJet->SetName(HistName);
        
        HistTempSingleHadronHard->Reset();
        HistTempSingleHadronSoft->Reset();
        //sprintf(HistName,"CountVspTSingleHadronSpectrumBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        //HistTempSingleHadron->SetName(HistName);
        
        fjInputs.resize(0);
        chargeList.resize(0);
        
        //temp hists for jets
        TH1D *HistTempJet2 = new TH1D("JetSpectrumBinTemp2", "Jet Spectrum pT", NpTJetBin, JetpTBin);
        TH1D *HistTempJet3 = new TH1D("JetSpectrumBinTemp3", "Jet Spectrum pT", NpTJetBin, JetpTBin);

        //temp hists for identified hadrons
        TH1D *tempPionsSoft = new TH1D("Soft Pion Spectrum Temp", "Pion Spectrum pT", NpTpionBin, piondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempPionsHard = new TH1D("Hard Pion Spectrum Temp", "Pion Spectrum pT", NpTpionBin, piondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempKaonsSoft = new TH1D("Soft Kaon Spectrum Temp", "Kaon Spectrum pT", NpTkaonBin, kaondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempKaonsHard = new TH1D("Hard Kaon Spectrum Temp", "Kaon Spectrum pT", NpTkaonBin, kaondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtonsSoft = new TH1D("Soft Proton Spectrum Temp", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtonsHard = new TH1D("Hard Proton Spectrum Temp", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtonsReco = new TH1D("Hard Proton Spectrum Temp", "Proton Spectrum pT", NpTprotonBin, protondata->GetXaxis()->GetXbins()->GetArray());

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        vector <fjcore::PseudoJet> UnSortedJets, SortedJets, UnSortedJets2, SortedJets2, UnSortedJets3, SortedJets3, constituents;
        
        //actually reading in, event loop
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
            Events++; //if(Events > 500) break;
            for(unsigned int i=0; i<hadrons.size(); i++){
                SN = i;
                PID= hadrons[i].get()->pid();
                E  = hadrons[i].get()->e();
                Px = hadrons[i].get()->px();
                Py = hadrons[i].get()->py();
                Pz = hadrons[i].get()->pz();
                Eta = hadrons[i].get()->eta();
                Y = hadrons[i].get()->rapidity();
                Phi = hadrons[i].get()->phi();
                pStat = hadrons[i].get()->pstat();
                mass = hadrons[i].get()->restmass();
                double PT = TMath::Sqrt( (Px*Px) + (Py*Py));
                
                if( fabs(Eta) < DetectorEtaCut && PT>0.01  &&  PID!=12 && PID!=14 && PID!=16 && PID!=18 ){
                    fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                    chargeList.push_back( pythia.particleData.charge( PID ) );
                }

                //strength for had hist filling
                float strength = 1.0f;
                
                //cutting for specific regimes
                if(k == 0 && PT > softend)
                    continue;
                if(k != 0 && PT < softend)
                    continue;
                if(k != 0 && PT > 1.1*stod(pTHatMax[k]))
                    continue;
                
                // Add this particle into SingleHadron spectrum
                if(fabs(Eta) < SingleHadronEtaCut && PT>0.01  && fabs(PID)>100 &&  pythia.particleData.charge(PID)!=0){
                    //cout<<PT<<" PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
                    if(PT < softend) HistTempSingleHadronSoft->Fill(PT,strength);
                    else HistTempSingleHadronHard->Fill(PT,strength);
                }
                
                //ID hadron spectra
                if(fabs(Y) < idHadronYCut){
                    if(abs(PID) == 211) {
                        if(PT < softend) tempPionsSoft->Fill(PT,strength);
                        else tempPionsHard->Fill(PT,strength);
                    }
                    if(abs(PID) == 321) {
                        if(PT < softend) tempKaonsSoft->Fill(PT,strength);
                        else tempKaonsHard->Fill(PT,strength);
                    }
                    if(abs(PID) == 2212) {
                        if(PT < softend) tempProtonsSoft->Fill(PT,strength);
                        else tempProtonsHard->Fill(PT,strength);
                        
                        if(pStat < 815) tempProtonsReco->Fill(PT,strength);
                    }
                } 
            }

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
                }
            }
            fjInputs.resize(0);
        }

        //xsec stuff
        double HardCrossSection = myfile->GetSigmaGen();
        double HardCrossSectionError = myfile->GetSigmaErr();
        if(k == 0) xsectotal = HardCrossSection; //set for first bin to match experimental value; end of reading cross section

        //cleaning up 0s in soft bin
        /*if(k == 0){
            for(int i=1; i <= tempPions->GetNbinsX(); i++)
                if(tempPions->GetBinContent(i) == 0 && tempPions->GetBinCenter(i) < softend) tempPions->SetBinContent(i,0.1);
            for(int i=1; i <= tempKaons->GetNbinsX(); i++)
                if(tempKaons->GetBinContent(i) == 0 && tempKaons->GetBinCenter(i) < softend) tempKaons->SetBinContent(i,0.1);
            for(int i=1; i <= tempProtons->GetNbinsX(); i++)
                if(tempProtons->GetBinContent(i) == 0 && tempProtons->GetBinCenter(i) < softend) tempProtons->SetBinContent(i,0.1);
        }*/

        //event count handling
        eventCount.push_back(Events);
        
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
        //for(int i = 0; i < NpTSingleHadronBin; i++)
            //if(HistTempSingleHadron->GetBinContent(i+1) < 5) HistTempSingleHadron->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTpionBin; i++)
            //if(tempPions->GetBinContent(i+1) < 2) tempPions->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTkaonBin; i++)
            //if(tempKaons->GetBinContent(i+1) < 2) tempKaons->SetBinContent(i+1,0);

        //for(int i = 0; i < NpTpionBin; i++)
            //if(tempProtons->GetBinContent(i+1) < 2) tempProtons->SetBinContent(i+1,0);
        
        //Write histogram into a root file
        HistTempJet->Sumw2(); HistTempJet->Write();
        HistTempJet2->Sumw2(); HistTempJet2->Write();
        HistTempJet3->Sumw2(); HistTempJet3->Write();
        tempPionsSoft->Sumw2(); tempPionsSoft->Write();
        tempPionsHard->Sumw2(); tempPionsHard->Write();
        tempKaonsSoft->Sumw2(); tempKaonsSoft->Write();
        tempKaonsHard->Sumw2(); tempKaonsHard->Write();
        tempProtonsSoft->Sumw2(); tempProtonsSoft->Write();
        tempProtonsHard->Sumw2(); tempProtonsHard->Write();
        tempProtonsReco->Sumw2(); tempProtonsReco->Write();
        HistTempSingleHadronSoft->Sumw2(); HistTempSingleHadronSoft->Write();
        HistTempSingleHadronHard->Sumw2(); HistTempSingleHadronHard->Write();
        
        //add to totals histograms
        double factor = HardCrossSection/(xsectotal*Events);
        HistTotalHadronSoft->Add(HistTempSingleHadronSoft,factor);
        HistTotalHadronHard->Add(HistTempSingleHadronHard,factor);
        HistTotalPionsSoft->Add(tempPionsSoft,factor);
        HistTotalPionsHard->Add(tempPionsHard,factor);
        HistTotalKaonsSoft->Add(tempKaonsSoft,factor);
        HistTotalKaonsHard->Add(tempKaonsHard,factor);
        HistTotalProtonsSoft->Add(tempProtonsSoft,factor);
        HistTotalProtonsHard->Add(tempProtonsHard,factor);
        HistRecoProtons->Add(tempProtonsReco,factor);
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

        //Save hadron pT hist as a png for convenience
        hadronComponents[k] = (TGraph*)GESingleHadron->Clone();
        hadronComponents[k]->SetLineColor(k+2);
        
        totalroot->cd();
    } //k-loop ends here (pTHatBin loop)
    
    delete HistTempJet; delete HistTempSingleHadronSoft; delete HistTempSingleHadronHard;

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

    //raw files
    HistTotalPionsSoft->Write("raw soft pions");
    HistTotalPionsHard->Write("raw hard pions");
    HistTotalKaonsSoft->Write("raw soft kaons");
    HistTotalKaonsHard->Write("raw hard kaons");
    HistTotalProtonsSoft->Write("raw soft protons");
    HistTotalProtonsHard->Write("raw hard protons");
    HistRecoProtons->Write("raw reco protons");

    //Scaling totals by global factors and the identified pions by bin centers
    HistTotalJet2->Scale(1.0/(2.0*JetEtaCut),"width");
    HistTotalJet3->Scale(1.0/(2.0*JetEtaCut),"width");
    scaleBins(HistTotalHadronSoft,(1.0/(2*M_PI*2.0*SingleHadronEtaCut)));
    scaleBins(HistTotalHadronHard,(1.0/(2*M_PI*2.0*SingleHadronEtaCut)));
    scaleBins(HistTotalPionsSoft,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistTotalPionsHard,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistTotalKaonsSoft,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistTotalKaonsHard,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistTotalProtonsSoft,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistTotalProtonsHard,(1.0/(2*M_PI*2.0*idHadronYCut)));
    scaleBins(HistRecoProtons,(1.0/(2*M_PI*2.0*idHadronYCut)));

    //hadrons
    //hadron data graph
    TMultiGraph *GEHadronTotal = new TMultiGraph();
    TFile hadron_file("/data/rjfgroup/rjf01/cameron.parker/data/HadronData.root");
    TDirectory* hadrondir = (TDirectory*)hadron_file.Get("Table 1");
    TGraphErrors* hadronData = (TGraphErrors*) hadrondir->Get("Graph1D_y1");
    TH1D* hadronDataHist = (TH1D*)hadrondir->Get("Hist1D_y1");
    hadronData->SetMarkerStyle(kCircle);
    hadronData->SetMarkerColor(kRed);
    hadronData->SetLineColor(kRed);
    hadronData->SetTitle("CMS");
    myRatioPlot(hadronData, HistTotalHadronSoft, HistTotalHadronHard, "Hadron Yields", true, true);

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
	hadronComp->Draw("alX");
    //hadronleg.DrawClone("Same");
	cTotalHadron->Print("plots/HadronTotalSpectrum.png");
    hadron_file.Close();
    
    HistFinalJet->Scale(1e6); //data listed in mb
    HistTotalJet2->Scale(1e6);
    HistTotalJet3->Scale(1e6);
	
    //create root file for total plots
    cout << "Creating ROOT output...";
    totalroot->cd();
    HistFinalJet->Write("jet radius 0.3");
    HistTotalHadronSoft->Write("soft hadrons");
    HistTotalHadronHard->Write("hard hadrons");
    TH1D* hadronhist = (TH1D*)HistTotalHadronSoft->Clone(); hadronhist->Add(HistTotalHadronHard); hadronhist->Write("hadrons");
    HistTotalJet2->Write("jet radius 0.2");
    HistTotalJet3->Write("jet radius 0.4");
    TH1D* pionhist = (TH1D*)HistTotalPionsSoft->Clone(); pionhist->Add(HistTotalPionsHard); pionhist->Write("smooth pions");
    HistTotalPionsSoft->Write("rough soft pions"); smoothBins(HistTotalPionsSoft); /*HistTotalPions->Smooth();*/ HistTotalPionsSoft->Write("smooth soft pions");
    HistTotalPionsHard->Write("rough hard pions"); smoothBins(HistTotalPionsHard); /*HistTotalPions->Smooth();*/ HistTotalPionsHard->Write("smooth hard pions");
    TH1D* kaonhist = (TH1D*)HistTotalKaonsSoft->Clone(); kaonhist->Add(HistTotalKaonsHard); kaonhist->Write("smooth kaons");
    HistTotalKaonsSoft->Write("rough soft kaons"); smoothBins(HistTotalKaonsSoft); /*HistTotalKaons->Smooth();*/ HistTotalKaonsSoft->Write("smooth soft kaons");
    HistTotalKaonsHard->Write("rough hard kaons"); smoothBins(HistTotalKaonsHard); /*HistTotalKaons->Smooth();*/ HistTotalKaonsHard->Write("smooth hard kaons");
    TH1D* protonhist = (TH1D*)HistTotalProtonsSoft->Clone(); protonhist->Add(HistTotalProtonsHard); protonhist->Write("smooth protons");
    HistTotalProtonsSoft->Write("rough soft protons"); smoothBins(HistTotalProtonsSoft); /*HistTotalProtons->Smooth();*/ HistTotalProtonsSoft->Write("smooth soft protons");
    HistTotalProtonsHard->Write("rough hard protons"); smoothBins(HistTotalProtonsHard); /*HistTotalProtons->Smooth();*/ HistTotalProtonsHard->Write("smooth hard protons");
    cout << "finished." << endl;

    //hadron graphs
    myRatioPlot(piongraph, HistTotalPionsSoft, HistTotalPionsHard, "Pion Yields", true, true);
    myRatioPlot(kaongraph, HistTotalKaonsSoft, HistTotalKaonsHard, "Kaon Yields", true, true);
    myRatioPlot(protongraph, HistTotalProtonsSoft, HistTotalProtonsHard, "Proton Yields", true, true);

    //Reco checks
    HistRecoProtons->Write("Reco proton spectrum");
    HistTotalProtonsSoft->Add(HistTotalProtonsHard);
    HistRecoProtons->Divide(HistTotalProtonsSoft);
    HistRecoProtons->Write("Reco proton ratio");

    //jets
    //jet data graph
    TMultiGraph *GEJetTotal = new TMultiGraph();
    TFile* jet_file = new TFile("/data/rjfgroup/rjf01/cameron.parker/data/JetData.root");
    TDirectory* jetdir = (TDirectory*)jet_file->Get("Table 4");
    TGraphErrors* jetData = (TGraphErrors*) jetdir->Get("Graph1D_y2");
    jetData->SetTitle("CMS");
    myRatioPlot(jetData,HistFinalJet,"Jet r3 Differential Cross sections",false,true);

    //r2 graph
    TGraphErrors* jetDataHist2 = (TGraphErrors*) jetdir->Get("Graph1D_y1");
    jetDataHist2->SetTitle("CMS");
    myRatioPlot(jetDataHist2,HistTotalJet2,"Jet r2 Differential Cross sections",false,true);

    //r4 graph
    TGraphErrors* jetDataHist4 = (TGraphErrors*) jetdir->Get("Graph1D_y3");
    jetDataHist4->SetTitle("CMS");
    myRatioPlot(jetDataHist4,HistTotalJet3,"Jet r4 Differential Cross sections",false,true);

    jet_file->Close();
    idhadron_file.Close();
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
