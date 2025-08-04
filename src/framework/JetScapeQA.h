/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef JETSCAPEQA_H
#define JETSCAPEQA_H

#include "JetScapeModuleBase.h"
#include <vector>
#include <Riostream.h>
#include <unordered_map>

#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

namespace Jetscape {

class JetScapeQA : public JetScapeModuleBase {
public:

    JetScapeQA() : JetScapeModuleBase() { SetId("JetScapeQA"); }
    JetScapeQA(string m_name) : JetScapeModuleBase(m_name) {SetId("JetScapeQA");}
    virtual ~JetScapeQA() {}
    
    virtual void Init();
    virtual void Exec();

    virtual void Finish();

private:

    TFile *fOutputFile = nullptr; // Output file for QA histograms
    int nEventsForQAHistograms = 100; // Number of events to fill
    bool enableEbyEQA = false; // Flag to enable/disable QA histograms
    string outputFileName; // Default output file name

    std::unordered_multimap<std::string,std::weak_ptr<JetScapeTask> > taskMap;
    void UpdateTaskMap();
    void PrintTasks();
    void PrintTaskMap();

    void DoEbyEQA();
    void DoQA();
    void WriteEbyEQA(string name, TH1 *h) {
        if (fOutputFile) {
            fOutputFile->cd();
            h->Write(name.c_str());
        }
    }

    void NormalizePerEvent();

    void JetPartonQA();
    void JetHadronQA();
    void JetPartonEbyEQA();
    //void JetHadronEbyEQA();

    //void HardProcessQA() {};
    void HardProcessEbyEQA();

    void SoftParticlizatonQA(); 
    //void SoftParticlizatonEbyEQA();

    void HydroEbyEQA() {};
    void ISEbyEQA() {};

    void PrintPDF();

    //Histograms ...
    TH1D *hJetPartonPt = nullptr; // Histogram for jet parton pT
    TH1D *hJetHadronPt = nullptr; // Histogram for jet hadron

    // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<JetScapeQA> reg;
};

} // namespace Jetscape
#endif
