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

#include "analysis.cc"

using namespace std;
using namespace Jetscape;

//using namespace Pythia8;
int main(int argc, char* argv[]){
    //Cut variables
    double DetectorEtaCut= 3.0;
    double pTmin = 5.0;
    int eventCount = 0;

    // For loop to open different pTHat bin files
    for (int k = 0; k<1; ++k){
        char HadronFile[300], pTBinString[100];
        sprintf(HadronFile,"/scratch/user/cameron.parker/projects/JETSCAPE/test/ko.dat.gz");
        auto myfile  = make_shared<JetScapeReaderAsciiGZ>(HadronFile);
     
        int  SN=0,PID=0;
        double Px, Py, Pz, E, Eta, Phi, pStat, mass, Y;
        int Events =0;
        int TriggeredJetNumber=0;

        //Data structures for events read in to save run time
        vector<shared_ptr<Hadron>> hadrons;
        
        //actually reading in
        while (!myfile->Finished()){
            myfile->Next();
            hadrons = myfile->GetHadrons();
            vector<shared_ptr<Hadron>> nucleons;
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
                    if(abs(PID) == 2212 || abs(PID) == 2112) nucleons.push_back(hadrons[i]);
                    if(PT > triggerPT){
                        triggerParticle = hadrons[i];
                        triggerPT = PT;
                    }
                } 
            }//End of hadron loop


        }
    } //k-loop ends here (pTHatBin loop)

    return 0;
}
