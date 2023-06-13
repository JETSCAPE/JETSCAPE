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
    if(stoi(getEcm()) != 5020){
        cout << "Center of mass energy does not match." << endl;
        return 0;
    }

    //get list of cross sections
    vector<vector<double>> xsecout = getXsecs(pTHatMin,pTHatMax);
	vector<double> xsecList = xsecout[0];
    vector<double> xsecErrorList = xsecout[1];
    double xsectotal = 0;
    for(int k = 0; k<NpTHardBin; k++) xsectotal += xsecList[k];
    //for(int k = 0; k<NpTHardBin; k++) cout << pTHatMin[k] << " " << pTHatMax[k] << " " << xsecList[k]*100000 << endl; //debugging line
    //cout << xsectotal << endl;

    //Cut variables
    double AssHadPtCut = 1.0;
    double deltaEtaCut = 1.0;
    double eYcut = 0.6;
    double triggerptcut[] = {4.0}; int nTrigRegions = 1;
    double assptmin[] = {1.0,2.0,5.0};  int nAssRegions = 3;
    double assptmax[] = {2.0,3.0,7.0};

    //output file for associated hadrons
    ofstream hadfile;
    hadfile.open("asshads.txt");
    
    //Variables for single hadron spectrum
    double ephiBin[17];
    for(int i = 0; i < 17; i++) ephiBin[i] = i*pi/16;
    double phibinw = ephiBin[1]-ephiBin[0];
    int NphieBin = sizeof(ephiBin)/sizeof(ephiBin[0])-1;
    TH1D *HistePhi[1][3]; //identified hadrons hists
    double totecount[1][3] = {0};
    string names[1][3];
    for(int i1 = 0; i1 < nTrigRegions; i1++){
        for(int i2 = 0; i2 < nAssRegions; i2++){
            names[i1][i2] = "D Azi Cor "+to_string(triggerptcut[i1])+"-"+to_string(triggerptcut[i1+1])+" "+
                to_string(assptmin[i2])+"-"+to_string(assptmax[i2]);
            HistePhi[i1][i2] = new TH1D("Hybrid Had. Prediction", names[i1][i2].c_str(), NphieBin, ephiBin);
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
        //sprintf(HadronFile,"test_out.dat");
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%s,%s) GeV",k,pTHatMin[k].c_str(),pTHatMax[k].c_str());
        
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int Events =0;
        int TriggeredJetNumber=0;
        int ecount[1][3] = {0};
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"root/SpectraBin%s_%s.root",pTHatMin[k].c_str(),pTHatMax[k].c_str());
        TFile* outFile = new TFile(outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];

        //temp hists for identified hadrons
        TH1D *epts = new TH1D("D counts", "D counts", 3, triggerptcut);
        TH1D *tempe[1][3]; //identified hadrons hists
        for(int i1 = 0; i1 < nTrigRegions; i1++)
            for(int i2 = 0; i2 < nAssRegions; i2++)
                tempe[i1][i2] = new TH1D(names[i1][i2].c_str(), names[i1][i2].c_str(), NphieBin, ephiBin);

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons, es, others;
        
        //actually reading in
        while (!myfile->Finished()){
            myfile->Next();
            hadrons = myfile->GetHadrons();

            //double counting protection
            if(k == 0){
                auto partons = myfile->GetPartonShowers();
                bool skip = false;

                for(int i = 0; i < partons.size(); i++){
                    if(partons[i]->GetPartonAt(0)->pt() > stod(pTHatMax[0])){
                        //cout << partons[i]->GetPartonAt(0)->pt() << endl;
                        skip = true;
                    }  
                }

                if(skip) continue;
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
                Eta = hadrons[i].get()->eta();
                Y = hadrons[i].get()->rapidity();
                Phi = hadrons[i].get()->phi();
                pStat = hadrons[i].get()->pstat();
                mass = hadrons[i].get()->restmass();
                double PT = TMath::Sqrt((Px*Px) + (Py*Py));                
                
                // making particle lists
                if(abs(PID) == 11){
                    if(abs(Y) < eYcut){
                        es.push_back(hadrons[i]);
                        epts->Fill(PT);
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
            for(int i = 0; i < es.size(); i++){
                for(int j = 0; j < others.size(); j++){
                    if(abs(es[i].get()->eta() - others[j].get()->eta()) < deltaEtaCut){
                        double deltaphi = abs(es[i].get()->phi() - others[j].get()->phi());
                        if(deltaphi > pi) deltaphi = (2*pi) - deltaphi;

                        for(int i1 = 0; i1 < nTrigRegions; i1++){
                            for(int i2 = 0; i2 < nAssRegions; i2++){
                                if(es[i].get()->pt() > triggerptcut[i1] && es[i].get()->pt() < triggerptcut[i1+1] 
                                    && others[j].get()->pt() > assptmin[i2] && others[j].get()->pt() < assptmax[i2]){
                                        tempe[i1][i2]->Fill(deltaphi, 1.0/(phibinw));
                                }
                            }
                        }
                    }
                }
            }

            //writing out associated hadrons
            if(es.size() > 0){
                hadfile << "pTHat: " << pTHatMin[k] << " to " << pTHatMax[k] << endl;
                for(int i = 0; i < others.size(); i++)
                    hadfile << others[i].get()->pid() << " " << others[i].get()->pt() << " " << others[i].get()->eta() << " " << others[i].get()->phi() << endl;
            }

            //clearing had vecs
            es.clear();
            others.clear();
        }
        
        //xsec and event count handling
        eventCount.push_back(Events);
        double HardCrossSection = xsecList[k];
        double HardCrossSectionError = xsecErrorList[k];
        
        //Dcounts
        for(int i1 = 0; i1 < nTrigRegions; i1++)
            for(int i2 = 0; i2 < nAssRegions; i2++)
                totecount[i1][i2] += epts->GetBinContent(i1+1)*HardCrossSection/(xsectotal*Events);
        
        //event info
        TVector EventInfo(3);
        EventInfo[0] = ecount[1][0];
        EventInfo[1] = ecount[1][1];
        EventInfo[2] = ecount[1][2];
        EventInfo.Write("EventInfo");
        
        //add to totals histograms
        //tempD->Scale(1.0/(Dcount));
        for(int i1 = 0; i1 < nTrigRegions; i1++){
            for(int i2 = 0; i2 < nAssRegions; i2++){
                tempe[i1][i2]->Write(names[i1][i2].c_str());
                HistePhi[i1][i2]->Add(tempe[i1][i2],HardCrossSection/(xsectotal*Events));
            }
        }

        myfile->Close();
        outFile->Close();
    } //k-loop ends here (pTHatBin loop)
    hadfile.close();

    //create root file for total plots
    TFile* totalroot = new TFile( "root/totals.root", "RECREATE");

    //Scaling totals by global factors and the identified pions by bin centers: dSigma/(2*pi*pT*dpT*dEta)
    //HistDPhi->SetNormFactor(1.0);
    for(int i1 = 0; i1 < nTrigRegions; i1++){
        for(int i2 = 0; i2 < nAssRegions; i2++){
            HistePhi[i1][i2]->Scale(1.0/totecount[i1][i2]);
            HistePhi[i1][i2]->GetXaxis()->SetTitle("Delta phi (rad)");
            HistePhi[i1][i2]->GetYaxis()->SetTitle("(1/Ne)(dNassc/dDelphi)");
 	        HistePhi[i1][i2]->Write(names[i1][i2].c_str());

            //plot output
            TCanvas *c = new TCanvas();
            HistePhi[i1][i2]->Draw();
            string plotname = "plots/" + names[i1][i2] + ".png";
            HistePhi[i1][i2]->GetYaxis()->SetRangeUser(0,10);
            c->SaveAs(plotname.c_str());
            c->Close();

            //ascii output
            ofstream asciifile;
            asciifile.open(names[i1][i2]);
            asciifile << "x\ty\tyerr" << endl;
            for(int i = 1; i <= HistePhi[i1][i2]->GetNbinsX(); i++)
                asciifile << HistePhi[i1][i2]->GetBinCenter(i) << "\t" << HistePhi[i1][i2]->GetBinContent(i)
                    << "\t" << HistePhi[i1][i2]->GetBinError(i) << endl;
            asciifile.close();

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
