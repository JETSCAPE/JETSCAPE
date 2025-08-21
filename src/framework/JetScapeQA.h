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
#include <set>
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

    // xml input options 
    TFile *fOutputFile = nullptr; // Output file for QA histograms
    int nEventByEventHistograms = 0;
    int nEventsForQAHistograms = 0;
    double th1_ptmax = -1; 
    double min_jetpt = -1.;
    int th1_nhadrons = -1;
    bool normalize_hgrams_per_event = false;
    string outputFileName; // Default output file name

    // keep a taskmap of all tasks in JetScape
    std::unordered_multimap<std::string,std::weak_ptr<JetScapeTask> > taskMap;
    void UpdateTaskMap();
    void PrintTasks();
    std::vector<std::tuple<int,string,string>> GetTaskInfo(bool sorted=true);
    void PrintTaskMap(bool sorted=true);

    // enumerator so that tasks selected for QA can
    // be routed to the proper type of QA (according
    // the base class)
    enum QA_TYPE {
        INITIAL_STATE, // <- trento
        HARD_PROCESS, // <- PythiaGun
        PREEQUILIBRIUM_DYNAMICS, // <- NullPreDynamics
        FLUID_DYNAMICS, // <- MUSIC
        JET_ENERGY_LOSS, // <- JetEnergyLoss
        HADRONIZATION,
        SOFT_PARTICLIZATION,
        NOT_IMPLEMENTED,
    };

    //based on input, select tasks for QA
    vector<std::tuple<string, QA_TYPE, vector<TH1*>>> qa_tasks;
    //based on input, select tasks for event-by-event QA
    vector<std::tuple<string, QA_TYPE, vector<TH1*>>> qa_EbyE_tasks;

    vector<TH1*> MakeHgrams(const string& task, QA_TYPE qa_type, int event=-1);
    void FillHgrams(std::tuple<string, QA_TYPE, vector<TH1*>>& qa);


    // flag because some histograms require information not available at 
    // JETSCAPE ::Init, and have to wait for the first JETSCAPE::Exec cycle
    bool has_first_exec = false;
    bool has_run_JEL = false; // Jet energy loss may have multiple tasks in list -- run all together

    void PrintPDF();

    static RegisterJetScapeModule<JetScapeQA> reg;
};
} // namespace Jetscape

#endif
