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
#include "TDirectory.h"

#include "analysis.cc"

using namespace std;
using namespace Jetscape;

//using namespace Pythia8;
int main(int argc, char* argv[]){
    //root initializations
    TApplication theApp("hist", &argc, argv);
    TFile* totalroot = new TFile( "/scratch/user/cameron.parker/projects/JETSCAPE/runs/ko/totals.root", "RECREATE");

    //Cut variables
    double DetectorEtaCut= 3.0;
    double dYcut= 0.5;
    double pTmin = 5.0;
    int eventCount = 0;

    //data calling
    TFile dataFile("/scratch/user/cameron.parker/projects/JETSCAPE/data/deuterons.root");
    TDirectory* trandir = (TDirectory*)dataFile.Get("Table 1 - Transverse");
    TH1D* trandata = (TH1D*) trandir->Get("Hist1D_y1");
    TGraphErrors* trangraph = (TGraphErrors*) trandir->Get("Graph1D_y1");
    TH1D *tranhist = getBlankCopy(trandata,"Transverse Spectrum","Transverse Spectrum");
    
    TDirectory* towarddir = (TDirectory*)dataFile.Get("Table 2 - Toward");
    TH1D* towarddata = (TH1D*) towarddir->Get("Hist1D_y1");
    TGraphErrors* towardgraph = (TGraphErrors*) towarddir->Get("Graph1D_y1");
    TH1D *towardhist = getBlankCopy(towarddata,"Transverse Spectrum","Transverse Spectrum");
    
    TDirectory* awaydir = (TDirectory*)dataFile.Get("Table 3 - Away");
    TH1D* awaydata = (TH1D*) awaydir->Get("Hist1D_y1");
    TGraphErrors* awaygraph = (TGraphErrors*) awaydir->Get("Graph1D_y1");
    TH1D *awayhist = getBlankCopy(awaydata,"Transverse Spectrum","Transverse Spectrum");

    // For loop to open different pTHat bin files
    int Events =0;
    for (int k = 0; k<1; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"/scratch/user/cameron.parker/projects/JETSCAPE/runs/ko/ko.dat.gz");
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
     
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int TriggeredJetNumber=0;

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        
        //actually reading in
        while (!myfile->Finished()){
            myfile->Next();
            hadrons = myfile->GetHadrons();
            vector<shared_ptr<Hadron>> nucleons, antinucleons;
            shared_ptr<Hadron> triggerParticle = hadrons[0];
            double triggerPT = 0.0;

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

                if(fabs(Y) < DetectorEtaCut){
                    //nucleons
                    if(PID == 2212 || PID == 2112){
                        nucleons.push_back(hadrons[i]);
                    }

                    //anti nucleons
                    if(PID == -2212 || PID == -2112){
                        antinucleons.push_back(hadrons[i]);
                    }

                    //triggers
                    if(PT > triggerPT){
                        triggerParticle = hadrons[i];
                        triggerPT = PT;
                    }
                } 
            }//End of hadron loop

            //if event kept
            if(triggerPT >= pTmin && (nucleons.size()>1 || antinucleons.size()>1)){
                eventCount++;
                
                for(int i=0; i<nucleons.size(); i++){
                    //skipping those out of bounds
                    if(fabs(nucleons[i].get()->rapidity()) < dYcut) continue;

                    //phi calculation
                    double deltaphi = fabs(triggerParticle.get()->phi() - nucleons[i].get()->phi());
                    if(deltaphi > pi) deltaphi = (2*pi) - deltaphi;

                    //hist filling on region
                    if(deltaphi < pi/3) towardhist->Fill(nucleons[i].get()->pt(), 0.5);
                    if(deltaphi < 2*pi/3 && deltaphi > pi/3) tranhist->Fill(nucleons[i].get()->pt(), 0.5);
                    if(deltaphi > 2*pi/3) awayhist->Fill(nucleons[i].get()->pt(), 0.5);
                }
            }
        }
    } //k-loop ends here (pTHatBin loop)

    //scaling
    double factor = 1.0/(Events*2.0*dYcut);
    tranhist->Scale(factor,"width");
    awayhist->Scale(factor,"width");
    towardhist->Scale(factor,"width");

    //end behavior
    totalroot->cd();
    tranhist->Write("Transverse");
    awayhist->Write("Away");
    towardhist->Write("Toward");
    totalroot->Close();
    cout << eventCount << " events kept" << endl;
    return 0;
}
