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
    double idHadronYCut = 1;
    double softend = 6.0;
    
    //reading data to get bins
    TFile dataroot( "/data/rjfgroup/rjf01/cameron.parker/data/LHC7000-ID-hads.root");
    TDirectory* piondir = (TDirectory*)dataroot.Get("Table 6");
    TDirectory* kaondir = (TDirectory*)dataroot.Get("Table 6");
    TDirectory* protondir = (TDirectory*)dataroot.Get("Table 6");

    //Variables for ID hadron hists
    TH1D *HistTotalPionsSoft = new TH1D("Pion Spectrum Soft", "Pion Spectrum pT", 22, 0.1, 1.2); //identified hadrons hists
    TH1D *HistTotalPionsHard = new TH1D("Pion Spectrum Hard", "Pion Spectrum pT", 22, 0.1, 1.2);
    TH1D *HistTotalKaonsSoft = new TH1D("Kaon Spectrum Soft", "Kaon Spectrum pT", 17, 0.2, 1.05);
    TH1D *HistTotalKaonsHard = new TH1D("Kaon Spectrum Hard", "Kaon Spectrum pT", 17, 0.2, 1.05);
    TH1D *HistTotalProtonsSoft = new TH1D("Proton Spectrum Soft", "Proton Spectrum pT", 27, 0.35, 1.7);
    TH1D *HistTotalProtonsHard = new TH1D("Proton Spectrum Hard", "Proton Spectrum pT", 27, 0.35, 1.7);
    cout << "Hadron graphs made" << endl;
    
    //doing the same for jets
    double JetpTBin[] = {20, 24, 28, 32, 38, 44, 50, 58, 66, 76, 86, 100}; //in GeV
    int NpTJetBin = sizeof(JetpTBin)/sizeof(JetpTBin[0])-1;     double JetpTMin = 10; //in GeV
    
    TFile jetdataroot( "/data/rjfgroup/rjf01/cameron.parker/data/LHC7000-jets.root");
    TDirectory* jetdir1 = (TDirectory*)jetdataroot.Get("Table 4");
    TH1D* jethist1 = new TH1D("R2 Jets", "R2 Jet Spectrum pT", NpTJetBin, JetpTBin);
    TDirectory* jetdir2 = (TDirectory*)jetdataroot.Get("Table 1");
    TH1D* jethist2 = new TH1D("R4 Jets", "R4 Jet Spectrum pT", NpTJetBin, JetpTBin);
    TDirectory* jetdir3 = (TDirectory*)jetdataroot.Get("Table 7");
    TH1D* jethist3 = new TH1D("R6 Jets", "R6 Jet Spectrum pT", NpTJetBin, JetpTBin);
    cout << "Jet graphs made" << endl;
    
    //jet def
    fjcore::JetDefinition jetDef1(fjcore::antikt_algorithm, 0.2);
    fjcore::JetDefinition jetDef2(fjcore::antikt_algorithm, 0.4);
    fjcore::JetDefinition jetDef3(fjcore::antikt_algorithm, 0.6);
    std::vector <fjcore::PseudoJet> SortedJets, UnsortedJets;

    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat.gz", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Y, Phi, pStat, mass, eta;
        int Events =0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"SpectraBin%s_%s",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        binfiles[k] = totalroot->mkdir(outFileName);
        binfiles[k]->cd();
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *tempPionsSoft = new TH1D("Soft Pion Spectrum Temp", "Pion Spectrum pT", 22, 0.1, 1.2);
        TH1D *tempPionsHard = new TH1D("Hard Pion Spectrum Temp", "Pion Spectrum pT", 22, 0.1, 1.2);
        TH1D *tempKaonsSoft = new TH1D("Soft Kaon Spectrum Temp", "Kaon Spectrum pT", 17, 0.2, 1.05);
        TH1D *tempKaonsHard = new TH1D("Hard Kaon Spectrum Temp", "Kaon Spectrum pT", 17, 0.2, 1.05);
        TH1D *tempProtonsSoft = new TH1D("Soft Proton Spectrum Temp", "Proton Spectrum pT", 27, 0.35, 1.7);
        TH1D *tempProtonsHard = new TH1D("Hard Proton Spectrum Temp", "Proton Spectrum pT", 27, 0.35, 1.7);
        TH1D* tempjethist1 = new TH1D("R2 Jets", "R2 Jet Spectrum pT", NpTJetBin, JetpTBin);
        TH1D* tempjethist2 = new TH1D("R4 Jets", "R4 Jet Spectrum pT", NpTJetBin, JetpTBin);
        TH1D* tempjethist3 = new TH1D("R6 Jets", "R6 Jet Spectrum pT", NpTJetBin, JetpTBin);

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

            //if(Events > 1000) break;

            //cout<<"Number of hadrons is: " << hadrons.size() << endl;
            Events++;
            std::vector <fjcore::PseudoJet> fjInputs;
            for(unsigned int i=0; i<hadrons.size(); i++){
                SN = i;
                PID= hadrons[i].get()->pid();
                E  = hadrons[i].get()->e();
                Px = hadrons[i].get()->px();
                Py = hadrons[i].get()->py();
                Pz = hadrons[i].get()->pz();
                eta = hadrons[i].get()->pseudorapidity();
                Y = hadrons[i].get()->rapidity();
                Phi = hadrons[i].get()->phi();
                pStat = hadrons[i].get()->pstat();
                mass = hadrons[i].get()->restmass();
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));
                
                if(PT>0.01 && PID!=12 && PID!=14 && PID!=16 && PID!=18 && PT > 0.15 && abs(eta) < 0.9){
                    fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                }      

                //cutting for specific regimes
                if(k == 0 && PT > softend)
                    continue;
                if(k != 0 && PT < softend)
                    continue;

                double strength = 1.0; //smoothing between smooth and hard transition          

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
                    }
                } 
            }

            //jet R=0.2
            fjcore::ClusterSequence clustSeq1(fjInputs, jetDef1);
            UnsortedJets = clustSeq1.inclusive_jets();
            SortedJets = sorted_by_pt(UnsortedJets);
            int pFast = SortedJets.size();
            for (auto jet: SortedJets){
                double jetpT = jet.perp();
                if(jetpT > stod(pTHatMax[k])*1.1) jetpT = stod(pTHatMax[k]); //catching high energy anomalies
                if(fabs(jet.pseudorapidity()) < 0.7) tempjethist1->Fill(jetpT);
            }

            //jet R=0.4
            fjcore::ClusterSequence clustSeq2(fjInputs, jetDef2);
            UnsortedJets = clustSeq2.inclusive_jets();
            SortedJets = sorted_by_pt(UnsortedJets);
            pFast = SortedJets.size();
            for (auto jet: SortedJets){
                double jetpT = jet.perp();
                if(jetpT > stod(pTHatMax[k])*1.1) jetpT = stod(pTHatMax[k]); //catching high energy anomalies
                if(fabs(jet.pseudorapidity()) < 0.5) tempjethist2->Fill(jetpT);
            }

            //jet R=0.6
            fjcore::ClusterSequence clustSeq3(fjInputs, jetDef3);
            UnsortedJets = clustSeq3.inclusive_jets();
            SortedJets = sorted_by_pt(UnsortedJets);
            pFast = SortedJets.size();
            for (auto jet: SortedJets){
                double jetpT = jet.perp();
                if(jetpT > stod(pTHatMax[k])*1.1) jetpT = stod(pTHatMax[k]); //catching high energy anomalies
                if(fabs(jet.pseudorapidity()) < 0.3) tempjethist3->Fill(jetpT);
            }
        }

        //xsec stuff
        double HardCrossSection = myfile->GetSigmaGen();
        double HardCrossSectionError =  myfile->GetSigmaErr();
        if(k == 0) xsectotal = HardCrossSection; //set for first bin to match experimental value; end of reading cross section
        
        //event count handling
        eventCount.push_back(Events);
        
        //Write histogram into a root file
        tempPionsSoft->Sumw2(); tempPionsSoft->Write();
        tempPionsHard->Sumw2(); tempPionsHard->Write();
        tempKaonsSoft->Sumw2(); tempKaonsSoft->Write();
        tempKaonsHard->Sumw2(); tempKaonsHard->Write();
        tempProtonsSoft->Sumw2(); tempProtonsSoft->Write();
        tempProtonsHard->Sumw2(); tempProtonsHard->Write();
        tempjethist1->Write();
        tempjethist2->Write();
        tempjethist3->Write();
        
        //add to totals histograms 
        HistTotalPionsSoft->Add(tempPionsSoft,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalPionsHard->Add(tempPionsHard,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalKaonsSoft->Add(tempKaonsSoft,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalKaonsHard->Add(tempKaonsHard,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalProtonsSoft->Add(tempProtonsSoft,HardCrossSection/(1.0*Events*xsectotal));
        HistTotalProtonsHard->Add(tempProtonsHard,HardCrossSection/(1.0*Events*xsectotal));
        jethist1->Add(tempjethist1,HardCrossSection*1.0e6/(1.0*Events));
        jethist2->Add(tempjethist2,HardCrossSection*1.0e6/(1.0*Events));
        jethist3->Add(tempjethist3,HardCrossSection*1.0e6/(1.0*Events));
		
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");

        totalroot->cd();
    } //k-loop ends here (pTHatBin loop)

    //raw files
    HistTotalPionsSoft->Write("raw soft pions");
    HistTotalPionsHard->Write("raw hard pions");
    HistTotalKaonsSoft->Write("raw soft kaons");
    HistTotalKaonsHard->Write("raw hard kaons");
    HistTotalProtonsSoft->Write("raw soft protons");
    HistTotalProtonsHard->Write("raw hard protons");

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    HistTotalPionsSoft->Scale(1./(2.0*idHadronYCut),"width");
    HistTotalPionsHard->Scale(1./(2.0*idHadronYCut),"width");
    HistTotalKaonsSoft->Scale(1./(2.0*idHadronYCut),"width");
    HistTotalKaonsHard->Scale(1./(2.0*idHadronYCut),"width");
    HistTotalProtonsSoft->Scale(1./(2.0*idHadronYCut),"width");
    HistTotalProtonsHard->Scale(1./(2.0*idHadronYCut),"width");
    jethist1->Scale(1.0/(2.0*0.7),"width"); //milli to nano
    jethist2->Scale(1.0/(2.0*0.5),"width");
    jethist3->Scale(1.0/(2.0*0.3),"width");
 	
    //create root file for total plots
    HistTotalPionsSoft->Write("rough soft pions"); smoothBins(HistTotalPionsSoft); /*HistTotalPions->Smooth();*/ HistTotalPionsSoft->Write("smooth soft pions");
    HistTotalPionsHard->Write("rough hard pions"); smoothBins(HistTotalPionsHard); /*HistTotalPions->Smooth();*/ HistTotalPionsHard->Write("smooth hard pions");
    TH1D* pionhist = (TH1D*)HistTotalPionsSoft->Clone(); pionhist->Add(HistTotalPionsHard); pionhist->Write("smooth pions");
    HistTotalKaonsSoft->Write("rough soft kaons"); smoothBins(HistTotalKaonsSoft); /*HistTotalKaons->Smooth();*/ HistTotalKaonsSoft->Write("smooth soft kaons");
    HistTotalKaonsHard->Write("rough hard kaons"); smoothBins(HistTotalKaonsHard); /*HistTotalKaons->Smooth();*/ HistTotalKaonsHard->Write("smooth hard kaons");
    TH1D* kaonhist = (TH1D*)HistTotalKaonsSoft->Clone(); kaonhist->Add(HistTotalKaonsHard); kaonhist->Write("smooth kaons");
    HistTotalProtonsSoft->Write("rough soft protons"); smoothBins(HistTotalProtonsSoft); /*HistTotalProtons->Smooth();*/ HistTotalProtonsSoft->Write("smooth soft protons");
    HistTotalProtonsHard->Write("rough hard protons"); smoothBins(HistTotalProtonsHard); /*HistTotalProtons->Smooth();*/ HistTotalProtonsHard->Write("smooth hard protons");
    TH1D* protonhist = (TH1D*)HistTotalProtonsSoft->Clone(); protonhist->Add(HistTotalProtonsHard); protonhist->Write("smooth protons");
    jethist1->Write("R = 0.2 jets"); smoothBins(jethist1); jethist1->Write("smooth R = 0.2 jets");
    jethist2->Write("R = 0.4 jets"); smoothBins(jethist2); jethist2->Write("smooth R = 0.4 jets");
    jethist3->Write("R = 0.6 jets"); smoothBins(jethist3); jethist3->Write("smooth R = 0.6 jets");
    totalroot->Close();

    //hadron graphs
    myRatioPlot((TGraphErrors*)piondir->Get("Graph1D_y1"), HistTotalPionsSoft, HistTotalPionsHard, "Pion Yields", true, true);
    myRatioPlot((TGraphErrors*)kaondir->Get("Graph1D_y2"), HistTotalKaonsSoft, HistTotalKaonsHard, "Kaon Yields", true, true);
    myRatioPlot((TGraphErrors*)protondir->Get("Graph1D_y3"), HistTotalProtonsSoft, HistTotalProtonsHard, "Proton Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir1->Get("Graph1D_y1"), jethist1, "R=0.2 Jet Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir2->Get("Graph1D_y1"), jethist2, "R=0.4 Jet Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir3->Get("Graph1D_y1"), jethist3, "R=0.6 Jet Yields", true, true);

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
