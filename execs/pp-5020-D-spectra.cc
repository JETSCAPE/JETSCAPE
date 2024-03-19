// lambda_c predictions for ALICE

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
#include <TH2.h>
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
    double recohads = 0;
    double fraghads = 0;

    //variable to total xsec
    double xsectotal = 0;

    //Cut variables
    double AssHadEtaCut = 0.8;
    double AssHadPtCut = 0.3;
    double deltaEtaCut = 1;
    double LYcut = 0.5;
    double triggerptcut[] = {0,3,5,8,16,24}; int ntrigbins = 5;
    double assptmin[] = {0.3,0.3,1,2,1}; int nassbins = 5;
    double assptmax[] = {10000,1,2,3,10000};
    
    //D spectra variables
    double spectraBins[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.0,3.35,3.8,4.4,5.1,6,7,8,9,10};
    int nSpectraBins = sizeof(spectraBins)/sizeof(spectraBins[0])-1;

    //Variables for single hadron spectrum
    double LphiBin[17];
    for(int i = 0; i < 17; i++) LphiBin[i] = i*pi/16;
    double phibinw = LphiBin[1]-LphiBin[0];
    int NphiLBin = sizeof(LphiBin)/sizeof(LphiBin[0])-1;
    TH1D *HistLPhi[ntrigbins][nassbins]; //identified hadrons hists
    double totLcount[ntrigbins] = {0};
    string names[ntrigbins][nassbins];
    for(int i1 = 0; i1 < ntrigbins; i1++){
        for(int i2 = 0; i2 < nassbins; i2++){
            names[i1][i2] = "D "+stringround(triggerptcut[i1],2)+"-"+stringround(triggerptcut[i1+1],2)+" GeV, hadrons "+
                stringround(assptmin[i2],2)+"-"+stringround(assptmax[i2],2)+" GeV";
            HistLPhi[i1][i2] = new TH1D("Hybrid Had. Prediction", names[i1][i2].c_str(), NphiLBin, LphiBin);
        }
    }   

    //graph declaration for adding hadron spectra
    TH1D *dRecoSpectra = new TH1D("Reco hads", "Reco hads;p_{T} (GeV);d^{2}#sigma/dp_{T}d#eta",nSpectraBins,spectraBins);
    TH1D *dFragSpectra = new TH1D("Frag hads", "Frag hads;p_{T} (GeV);d^{2}#sigma/dp_{T}d#eta",nSpectraBins,spectraBins);
    TH1D *HistTotalHadron = new TH1D("HadronSpectrumBin", "Combined Hadron pT Spectrum;p_{T} (GeV);d^{2}#sigma/dp_{T}d#eta", nSpectraBins, spectraBins); //Total hist for hadrons

    //debugging
    bool flag = true;
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    int Events =0;
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int TriggeredJetNumber=0;
        int Lcount[ntrigbins][nassbins] = {0};
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s.root",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        TFile* outFile = new TFile(outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *Lpts = new TH1D("D counts", "D counts", ntrigbins, triggerptcut);
        TH1D *tempL[ntrigbins][nassbins]; //identified hadrons hists
        for(int i1 = 0; i1 < ntrigbins; i1++)
            for(int i2 = 0; i2 < nassbins; i2++)
                tempL[i1][i2] = new TH1D(names[i1][i2].c_str(), names[i1][i2].c_str(), NphiLBin, LphiBin);

        //hists for reco vs fragmentation
        TH2D *recoHist[ntrigbins][nassbins]; //identified hadrons hists
        for(int i1 = 0; i1 < ntrigbins; i1++)
            for(int i2 = 0; i2 < nassbins; i2++)
                recoHist[i1][i2] = new TH2D(names[i1][i2].c_str(), names[i1][i2].c_str(), 2, 10, 30, 2, 10, 30);

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons, Ls, others;
        
        //actually reading in
        while (!myfile->Finished()){
            myfile->Next();
            //cout << "Geting hadrons" << endl;
            hadrons = myfile->GetHadrons();

            //cout<<"Number of hadrons is: " << hadrons.size() << endl;
            Events++;
            //if(Events == 10000) break;

            //particle loop
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
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));       

                if(pStat == 811) recohads++;
                else fraghads++;         
                
                // making particle lists
                if(abs(PID) == 411 || abs(PID) == 413 || abs(PID) == 421){
                    if(pStat==811) dRecoSpectra->Fill(PT);
                    if(pStat==821) dFragSpectra->Fill(PT);
                    if(abs(Y) < LYcut){
                        Ls.push_back(hadrons[i]);
                        Lpts->Fill(PT);
                    }
                }
                
                //associated particles are electrons, muons, pions, kaons, protons
                if(abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 211 || abs(PID) == 321 || abs(PID) == 2212){
                    if(PT > AssHadPtCut){
                        others.push_back(hadrons[i]);
                    }

                    HistTotalHadron->Fill(PT);
                }
                
            }

            //event variables

            //forming and binning pairs
            for(int i = 0; i < Ls.size(); i++){
                for(int j = 0; j < others.size(); j++){
                    if(abs(Ls[i].get()->eta() - others[j].get()->eta()) < deltaEtaCut){
                        double deltaphi = abs(Ls[i].get()->phi() - others[j].get()->phi());
                        if(deltaphi > pi) deltaphi = (2*pi) - deltaphi;

                        for(int i1 = 0; i1 < ntrigbins; i1++){
                            for(int i2 = 0; i2 < nassbins; i2++){
                                if(Ls[i].get()->pt() > triggerptcut[i1] && Ls[i].get()->pt() < triggerptcut[i1+1] 
                                    && others[j].get()->pt() > assptmin[i2] && others[j].get()->pt() < assptmax[i2]){
                                        tempL[i1][i2]->Fill(deltaphi, 1.0/(phibinw));
                                        recoHist[i1][i2]->Fill(Ls[i].get()->pstat()-800, others[j].get()->pstat()-800);
                                        //cout << Ls[i].get()->pstat() << " " << others[j].get()->pstat() << endl;
                                }
                            }
                        }
                    }
                }
            }

            //clearing had vecs
            if(Ls.size() == 0) flag = false;
            Ls.clear();
            others.clear();
            //cout << " . Event over." << endl;
        }
        
        //xsec and event count handling
        eventCount.push_back(Events);
        double HardCrossSection = myfile->GetSigmaGen();
        //cout << "Xsec: " << HardCrossSection << endl;
        double HardCrossSectionError = myfile->GetSigmaErr();
        xsectotal += HardCrossSection;
        
        //Dcounts
        for(int i1 = 0; i1 < ntrigbins; i1++)
            totLcount[i1] += Lpts->GetBinContent(i1+1);
        
        //event info
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");
        
        //add to totals histograms
        //tempD->Scale(1.0/(Dcount));
        for(int i1 = 0; i1 < ntrigbins; i1++){
            for(int i2 = 0; i2 < nassbins; i2++){
                tempL[i1][i2]->GetXaxis()->SetTitle("D Reco Label");
                tempL[i1][i2]->GetYaxis()->SetTitle("Hadron Reco Label");
                tempL[i1][i2]->Write(names[i1][i2].c_str());
                recoHist[i1][i2]->Write(names[i1][i2].c_str());
                HistLPhi[i1][i2]->Add(tempL[i1][i2]);
            }
        }

        myfile->Close();
        outFile->Close();
    } //k-loop ends here (pTHatBin loop)

    //Final writing
    TFile hadron_file("/scratch/user/cameron.parker/projects/JETSCAPE/data/D-meson-data.root");
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");
    cout << "Got data file" << endl;

    THStack *histStack = new THStack("Reco vs. Frag", "Reco vs. Frag");
    scaleBins(dRecoSpectra,1.0/Events); histStack->Add(dRecoSpectra);
    scaleBins(dFragSpectra,1.0/Events); histStack->Add(dFragSpectra);
    dRecoSpectra->Write("Reco Spectrum");
    dFragSpectra->Write("Frag Spectrum");
    histStack->Write("Revo vs Frag");

    scaleBins(HistTotalHadron,1.0/Events);
    HistTotalHadron->Write("pT Spectrum");

    int i = 1;
    for(int i2 = 0; i2 < nassbins; i2++){
        for(int i1 = 0; i1 < ntrigbins; i1++){
            HistLPhi[i1][i2]->Scale(1.0/(2*totLcount[i1]));
            HistLPhi[i1][i2]->GetXaxis()->SetTitle("Delta phi (rad)");
            HistLPhi[i1][i2]->GetYaxis()->SetTitle("(1/N_D)(dN_{assc}/d#Delta#phi)");
 	        HistLPhi[i1][i2]->Write(names[i1][i2].c_str());

            //data comparison
            if(i1==0 || i2==4) continue;
            string tablename = "Table " + to_string(i); i++;
            TDirectory* hadrondir = (TDirectory*)hadron_file.Get(tablename.c_str());
            TGraphErrors* hadronData = (TGraphErrors*) hadrondir->Get("Graph1D_y1");
            TH1D* temphist = getZeroedHist(HistLPhi[i1][i2]);
            temphist->GetXaxis()->SetRangeUser(0,pi);
            string xname = "#Delta#phi", yname = "1/N_{D} dN_{pairs}/#Delta#phi";
            ratioPlot(hadronData,temphist,names[i1][i2],xname,yname);
        }
    }
    totalroot->Close();

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    cout<<"Recombined hadrons: "<<recohads<<endl;
    cout<<"Fragmented hadrons: "<<fraghads<<endl;
    
    //test comment
    //debugging
    //cout << flag << endl;

    return 0;
}
