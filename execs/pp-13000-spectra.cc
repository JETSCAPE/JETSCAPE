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
    TH1D *HistTotalPionsSoft = new TH1D("Pion Spectrum Soft", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray()); //identified hadrons hists
    TH1D *HistTotalPionsHard = new TH1D("Pion Spectrum Hard", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalKaonsSoft = new TH1D("Kaon Spectrum Soft", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalKaonsHard = new TH1D("Kaon Spectrum Hard", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalProtonsSoft = new TH1D("Proton Spectrum Soft", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());
    TH1D *HistTotalProtonsHard = new TH1D("Proton Spectrum Hard", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());
    
    //doing the same for jets
    TFile jetdataroot( "/data/rjfgroup/rjf01/cameron.parker/data/LHC13000-jets.root");
    TDirectory* jetdir1 = (TDirectory*)jetdataroot.Get("Table 1"); TH1D* jetdata1 = (TH1D*)jetdir1->Get("Hist1D_y1");
    TH1D* jethist1 = getBlankCopy(jetdata1,"Jets low y","Jets low y");
    TDirectory* jetdir2 = (TDirectory*)jetdataroot.Get("Table 2"); TH1D* jetdata2 = (TH1D*)jetdir2->Get("Hist1D_y1");
    TH1D* jethist2 = getBlankCopy(jetdata2,"Jets mid y","Jets mid y");
    TDirectory* jetdir3 = (TDirectory*)jetdataroot.Get("Table 3"); TH1D* jetdata3 = (TH1D*)jetdir3->Get("Hist1D_y1");
    TH1D* jethist3 = getBlankCopy(jetdata3,"Jets high y","Jets high y");
    
    //jet def
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, 0.4);
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
        TH1D *tempPionsSoft = new TH1D("Soft Pion Spectrum Temp", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempPionsHard = new TH1D("Hard Pion Spectrum Temp", "Pion Spectrum pT", piondata->GetNbinsX(), piondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempKaonsSoft = new TH1D("Soft Kaon Spectrum Temp", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempKaonsHard = new TH1D("Hard Kaon Spectrum Temp", "Kaon Spectrum pT", kaondata->GetNbinsX(), kaondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtonsSoft = new TH1D("Soft Proton Spectrum Temp", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());
        TH1D *tempProtonsHard = new TH1D("Hard Proton Spectrum Temp", "Proton Spectrum pT", protondata->GetNbinsX(), protondata->GetXaxis()->GetXbins()->GetArray());
        TH1D* tempjethist1 = getBlankCopy(jetdata1,"Temp Jets low y","Temp Jets low y");
        TH1D* tempjethist2 = getBlankCopy(jetdata2,"Temp Jets mid y","Temp Jets mid y");
        TH1D* tempjethist3 = getBlankCopy(jetdata3,"Temp Jets high y","Temp Jets high y");

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
            std::vector <fjcore::PseudoJet> fjInputs;
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
                
                if(PT>0.01 && PID!=12 && PID!=14 && PID!=16 && PID!=18){
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

            //jet calcs
            fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
            UnsortedJets = clustSeq.inclusive_jets(100.);
            SortedJets = sorted_by_pt(UnsortedJets);
            int pFast = SortedJets.size();
            for (auto jet: SortedJets){
                double jetpT = jet.perp();
                if(jetpT > stod(pTHatMax[k])*1.1) jetpT = stod(pTHatMax[k]); //catching high energy anomalies
                if(fabs(jet.rapidity()) < 0.5) tempjethist1->Fill(jetpT);
                if(fabs(jet.rapidity()) > 0.5 and fabs(jet.rapidity()) < 1.0) tempjethist2->Fill(jetpT);
                if(fabs(jet.rapidity()) > 1.0 and fabs(jet.rapidity()) < 1.5) tempjethist3->Fill(jetpT);
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
        jethist1->Add(tempjethist1,HardCrossSection/(1.0*Events));
        jethist2->Add(tempjethist2,HardCrossSection/(1.0*Events));
        jethist3->Add(tempjethist3,HardCrossSection/(1.0*Events));
		
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
    jethist1->Scale(2000000000.0,"width"); //milli to pico times 2 for half the rapidity range
    jethist2->Scale(2000000000.0,"width");
    jethist3->Scale(2000000000.0,"width");
 	
    //create root file for total plots
    HistTotalPionsSoft->Write("rough soft pions"); smoothBins(HistTotalPionsSoft); /*HistTotalPions->Smooth();*/ HistTotalPionsSoft->Write("smooth soft pions");
    HistTotalPionsHard->Write("rough hard pions"); smoothBins(HistTotalPionsHard); /*HistTotalPions->Smooth();*/ HistTotalPionsHard->Write("smooth hard pions");
    HistTotalKaonsSoft->Write("rough soft kaons"); smoothBins(HistTotalKaonsSoft); /*HistTotalKaons->Smooth();*/ HistTotalKaonsSoft->Write("smooth soft kaons");
    HistTotalKaonsHard->Write("rough hard kaons"); smoothBins(HistTotalKaonsHard); /*HistTotalKaons->Smooth();*/ HistTotalKaonsHard->Write("smooth hard kaons");
    HistTotalProtonsSoft->Write("rough soft protons"); smoothBins(HistTotalProtonsSoft); /*HistTotalProtons->Smooth();*/ HistTotalProtonsSoft->Write("smooth soft protons");
    HistTotalProtonsHard->Write("rough hard protons"); smoothBins(HistTotalProtonsHard); /*HistTotalProtons->Smooth();*/ HistTotalProtonsHard->Write("smooth hard protons");
    jethist1->Write("low y jets"); smoothBins(jethist1); jethist1->Write("smooth low y jets");
    jethist2->Write("mid y jets"); smoothBins(jethist2); jethist2->Write("smooth mid y jets");
    jethist3->Write("high y jets"); smoothBins(jethist3); jethist3->Write("smooth high y jets");
    totalroot->Close();

    //hadron graphs
    myRatioPlot((TGraphErrors*)piondir->Get("Graph1D_y1"), HistTotalPionsSoft, HistTotalPionsHard, "Pion Yields", true, true);
    myRatioPlot((TGraphErrors*)kaondir->Get("Graph1D_y1"), HistTotalKaonsSoft, HistTotalKaonsHard, "Kaon Yields", true, true);
    myRatioPlot((TGraphErrors*)protondir->Get("Graph1D_y1"), HistTotalProtonsSoft, HistTotalProtonsHard, "Proton Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir1->Get("Graph1D_y1"), jethist1, "Low y Jet Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir2->Get("Graph1D_y1"), jethist2, "Mid y Jet Yields", true, true);
    myRatioPlot((TGraphErrors*)jetdir3->Get("Graph1D_y1"), jethist3, "High y Jet Yields", true, true);

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
