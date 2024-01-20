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
#include <TFile.h>
#include <TVector.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TH2D.h"

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

    //get list of cross sections
    double xsectotal = 0;

    //Cut variables
    double AssHadEtaCut = 0.8;
    double AssHadPtCut = 0.3;
    double deltaEtaCut = 1;
    double LYcut = 0.5;
    double triggerptcut[] = {3,5,8,16}; int ntrigbins = 3;
    double assptmin[] = {0.3,0.3,1}; int nassbins = 3;
    double assptmax[] = {10000,1,10000};
    
    //Variables for single hadron spectrum
    double LphiBin[17];
    for(int i = 0; i < 17; i++) LphiBin[i] = i*pi/16;
    double phibinw = LphiBin[1]-LphiBin[0];
    int NphiLBin = sizeof(LphiBin)/sizeof(LphiBin[0])-1;
    TH1D *tempL[3][3]; //identified hadrons hists
    double totLcount[3][3] = {0};
    string names[3][3];
    for(int i1 = 0; i1 < 3; i1++){
        for(int i2 = 0; i2 < 3; i2++){
            names[i1][i2] = "L Azi Cor "+to_string(triggerptcut[i1])+"-"+to_string(triggerptcut[i1+1])+" "+
                to_string(assptmin[i2])+"-"+to_string(assptmax[i2]);
            tempL[i1][i2] = new TH1D("Hybrid Had. prediction", names[i1][i2].c_str(), NphiLBin, LphiBin);
        }
    }   

    //graph declaration for adding hadron spectra
    TMultiGraph* hadronComp = new TMultiGraph();
    TGraph* hadronComponents[NpTHardBin];
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    for (int k = 0; k<NpTHardBin; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat", pTHatMin[k].c_str(), pTHatMax[k].c_str());
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int Events =0;
        int TriggeredJetNumber=0;
        int Lcount[3][3] = {0};
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s.root",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        TFile* outFile = new TFile(outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *Lpts = new TH1D("L counts", "L counts", 3, triggerptcut);
        
        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons, Ls, others;

        //hists for reco vs fragmentation
        TH2D *recoHist[ntrigbins][nassbins]; //identified hadrons hists
        for(int i1 = 0; i1 < ntrigbins; i1++)
            for(int i2 = 0; i2 < nassbins; i2++)
                recoHist[i1][i2] = new TH2D(names[i1][i2].c_str(), names[i1][i2].c_str(), 2, 10, 30, 2, 10, 30);
        
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
                Y = hadrons[i].get()->rapidity();
                Phi = hadrons[i].get()->phi();
                pStat = hadrons[i].get()->pstat();
                mass = hadrons[i].get()->restmass();
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));                
                
                // making particle lists
                if(abs(PID) == 4122){
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
                }
                
            }

            //forming and binning pairs
            for(int i = 0; i < Ls.size(); i++){
                for(int j = 0; j < others.size(); j++){
                    if(abs(Ls[i].get()->eta() - others[j].get()->eta()) < deltaEtaCut){
                        double deltaphi = abs(Ls[i].get()->phi() - others[j].get()->phi());
                        if(deltaphi > pi) deltaphi = (2*pi) - deltaphi;

                        for(int i1 = 0; i1 < 3; i1++){
                            for(int i2 = 0; i2 < 3; i2++){
                                if(Ls[i].get()->pt() > triggerptcut[i1] && Ls[i].get()->pt() < triggerptcut[i1+1] 
                                    && others[j].get()->pt() > assptmin[i2] && others[j].get()->pt() < assptmax[i2]){
                                        tempL[i1][i2]->Fill(deltaphi);
                                        recoHist[i1][i2]->Fill(Ls[i].get()->pstat()-800, others[j].get()->pstat()-800);
                                }
                            }
                        }
                    }
                }
            }

            //clearing had vecs
            Ls.clear();
            others.clear();
        }
        
        //xsec and event count handling
        eventCount.push_back(Events);
        double HardCrossSection = myfile->GetSigmaGen();
        //cout << "Xsec: " << HardCrossSection << endl;
        double HardCrossSectionError = myfile->GetSigmaErr();
        xsectotal += HardCrossSection;
        
        //Dcounts
        for(int i1 = 0; i1 < 3; i1++)
            for(int i2 = 0; i2 < 3; i2++)
                totLcount[i1][i2] += Lpts->GetBinContent(i1+1);
        
        //event info
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");
        
        //add to totals histograms
        //tempD->Scale(1.0/(Dcount));
        for(int i1 = 0; i1 < 3; i1++){
            for(int i2 = 0; i2 < 3; i2++){
                tempL[i1][i2]->Write(names[i1][i2].c_str());
                recoHist[i1][i2]->Write(names[i1][i2].c_str());
            }
        }

        myfile->Close();
        outFile->Close();
    } //k-loop ends here (pTHatBin loop)

    //create root file for total plots
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    //HistDPhi->SetNormFactor(1.0);
    for(int i1 = 0; i1 < 3; i1++){
        for(int i2 = 0; i2 < 3; i2++){
            tempL[i1][i2]->Scale(1.0/(2*totLcount[i1][i2]),"width");
            tempL[i1][i2]->GetXaxis()->SetTitle("#Delta#phi (rad)");
            tempL[i1][i2]->GetYaxis()->SetTitle("(1/N_{#Lambda} dN_{pairs}/d#Delta#phi)");
 	        tempL[i1][i2]->Write(names[i1][i2].c_str());

            //plot output
            TCanvas *c = new TCanvas();
            tempL[i1][i2]->Draw();
            string plotname = "plots/" + names[i1][i2] + ".png";
            //tempL[i1][i2]->GetYaxis()->SetRangeUser(0,10);
            c->SaveAs(plotname.c_str());
            c->Close();
        }
    }

    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Program run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    
    //test comment
    //debugging

    return 0;
}
