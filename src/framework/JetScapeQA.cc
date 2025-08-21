#include "JetScapeQA.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include <sstream>
#include "tinyxml2.h"
#include "JetScapeSignalManager.h"
#include "TrentoInitial.h"
#include "iSpectraSamplerWrapper.h"
#include "PythiaGun.h"

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<JetScapeQA> JetScapeQA::reg("JetScapeQA");

void JetScapeQA::Init() {
    JetScapeModuleBase::Init();
    JSINFO << "Initialize JetScapeQA : " << GetId() << " ...";

    // eveny-by-event qa is lazy -- it will only do the types selected
    // user can enter "none" or "all"
    auto ebye_option = GetXMLElementText({"JetScapeQA", "EventByEventTaskList"});
    bool ebye_all  = false;
    bool ebye_none  = false;
    vector<string> tasks_EbyE {};

    std::istringstream iss(ebye_option);
    string word;
    while (iss >> word) {
      if (word == "none") {
        ebye_none = true;
        break;
      } else if( word == "all") {
        ebye_all = true;
        break;
      } else {
        tasks_EbyE.push_back(word);
        JSINFO << " Doing event-by-event QA for task: " << word;
      }
    }
    // }

    // qa (not the event-by-event ones) are greedy -- they will do all
    // tasks except those excempted
    // user can enter "none" or "all"
    auto skip_qa_option = GetXMLElementText({"JetScapeQA", "SkipQAList"});
    bool skip_qa_all = false;
    bool skip_qa_none = false;
    vector<string> skip_qa_tasks{};
    iss = std::istringstream(skip_qa_option);
    while (iss >> word) {
      if (word == "all") {
        skip_qa_all = true;
        break;
      } else if (word == "none") {
        skip_qa_none = true;
        break;
      } else {
        skip_qa_tasks.push_back(word);
        JSINFO << " Skipping implemented QA module for task: " << word;
      }
    }

    // read in number of events used in event-by-event and in qa
    auto outputFileName = GetXMLElementText({"JetScapeQA", "outputFileName"}, true);
    JSINFO << "QA output FileName: " << outputFileName;
    nEventByEventHistograms = GetXMLElementInt({"JetScapeQA", "nEventByEventHistograms"});
    nEventsForQAHistograms = GetXMLElementInt({"JetScapeQA", "nEventsForQAHistograms"});
    normalize_hgrams_per_event = (GetXMLElementInt({"JetScapeQA", "normalizeHgramsPerEvent"}, false) > 0);

    // upper-bound for ptmax of histograms
    th1_ptmax = GetXMLElementDouble({"JetScapeQA","th1_ptmax"}, true);
    min_jetpt = GetXMLElementDouble({"JetScapeQA","min_jetpt"}, true);
    th1_nhadrons = GetXMLElementInt({"JetScapeQA","th1_nhadrons"}, true);

    JSINFO << "JetScapeQA: th1_ptmax = " << th1_ptmax;
    JSINFO << "JetScapeQA: number of event-by-event histograms = " << nEventByEventHistograms;
    JSINFO << "JetScapeQA: number of events used in QA histograms = " << nEventsForQAHistograms;

    fOutputFile=new TFile(outputFileName.c_str(), "RECREATE");

    // fill the task lists of QA to do
    UpdateTaskMap();
    auto taskInfo = GetTaskInfo();
    for (auto& task_tuple : taskInfo) {
      const auto task = std::get<1>(task_tuple);
      JSINFO << " DEBUG " << std::get<0>(task_tuple) << " " << std::get<1>(task_tuple) << " || " << std::get<2>(task_tuple);

      // find the base class of the task
      auto qa_type = QA_TYPE::NOT_IMPLEMENTED;
      if (task == "Trento") {
          qa_type = QA_TYPE::INITIAL_STATE;
      } else if (task == "PythiaGun") {
          qa_type = QA_TYPE::HARD_PROCESS;
      } else if (task == "Hadronization") {
          qa_type = QA_TYPE::HADRONIZATION;
      } else if (task == "NullPreDynamics") {
          qa_type = QA_TYPE::PREEQUILIBRIUM_DYNAMICS;
      } else if (task == "MUSIC") {
            qa_type = QA_TYPE::FLUID_DYNAMICS;
      } else if (task == "iSS") {
            qa_type = QA_TYPE::SOFT_PARTICLIZATION;
      } else if (task == "JetEnergyLoss") {
            qa_type = QA_TYPE::JET_ENERGY_LOSS;
      }
        JSINFO << " TRYING TO ADD TASK " << task << " to EbyE QA";
        JSINFO << " ALSO " << skip_qa_none ;

      if (qa_type == QA_TYPE::NOT_IMPLEMENTED) { continue; }

      // add it to EbyE tasks
      if (!ebye_none) {
        if (ebye_all || std::find(tasks_EbyE.begin(), tasks_EbyE.end(), task) != tasks_EbyE.end()) {
          JSINFO << " Adding task " << task << " to EbyE QA";
          qa_EbyE_tasks.push_back({task, qa_type,{}}); // make new histograms every event
        }
      }

      if (skip_qa_all) continue;
      bool noskip = std::find(skip_qa_tasks.begin(), skip_qa_tasks.end(), task) == skip_qa_tasks.end();
      if (skip_qa_none || noskip) {
        // Add the task to the QA list
        JSINFO << " DEBUG :: Adding task " << task;
        qa_tasks.push_back({task, qa_type, MakeHgrams(task,qa_type)});
      }
    }

    PrintTaskMap();
    JSINFO << "   -- listings of JetScapeQA tasks selected for QA histograms --";
    JSINFO << " JetScapeQA Tasks: (name, n-histograms at initialization)";
    for (auto& qa : qa_tasks) {
      JSINFO << "   + " <<  std::left << std::setw(20) <<  std::get<0>(qa) << " " <<std::get<2>(qa).size();
    }
    JSINFO << " Event-by-event JetScapeQA Tasks: (name, n-hgrams at initialization)";
    for (auto& qa : qa_EbyE_tasks) {
      JSINFO << "   + " <<  std::left << std::setw(20) <<  std::get<0>(qa) << " " <<std::get<2>(qa).size();
    }

    // initialize the histograms for the different kinds of input
}


void JetScapeQA::Exec() {
  has_first_exec = true;
  has_run_JEL = false;
  VERBOSE(8) << "Executing JetScapeQA : " << GetId() << " ...";
  UpdateTaskMap();
  //PrintTasks();
  //PrintTaskMap();
  int current_event = GetCurrentEvent();

  if (current_event < nEventByEventHistograms) {
    for (auto &qa : qa_EbyE_tasks) {
      std::get<2>(qa) = MakeHgrams(std::get<0>(qa), std::get<1>(qa), current_event);
      FillHgrams(qa);
      for (auto &h : std::get<2>(qa)) {
        if (fOutputFile) {
          // fOutputFile->cd();
          h->Write();
        }
        delete h;
      }
      std::get<2>(qa).clear();
    }
  }

  if (current_event < nEventsForQAHistograms) {
    for (auto &qa : qa_tasks) {
      FillHgrams(qa);
    }
  }

}

vector<TH1*> JetScapeQA::MakeHgrams(const string& task, QA_TYPE qa_type, int event) {
    vector<TH1*> hgrams;

    // see if there is anthing to get

    auto it = taskMap.find(task);
    string etag = (event >= 0) ? Form("_event_%d", event) : "";

    switch (qa_type) {
    case QA_TYPE::INITIAL_STATE: {
      auto istate = std::dynamic_pointer_cast<InitialState>(it->second.lock());
      float xstep = istate->GetXStep();
      float xmax = istate->GetXMax();
      float ystep = istate->GetYStep();
      float ymax = istate->GetYMax();
      int nx = 2*xmax / xstep;
      int ny = 2*ymax / ystep;

      // get a 2D histogram for entropy
      hgrams.push_back(new TH2D(Form("%s_entropy%s", task.c_str(),etag.c_str()), Form("%s Entropy;X [fm];Y [fm]", task.c_str()), nx, -xmax, xmax, ny, -ymax, ymax));
      hgrams.push_back(new TH1D(Form("%s_centrality%s", task.c_str(),etag.c_str()), Form("%s centrality; centrality [%%];", task.c_str()), 100, 0., 100.));
      // add additional histograms as wanted; see class InitialState.h for options
    }
    break;

    case QA_TYPE::HARD_PROCESS: {
      auto hp = std::dynamic_pointer_cast<HardProcess>(it->second.lock());
      auto pthat = hp->GetPtHat();
      hgrams.push_back(new TH1D( Form("%s_parton_pt%s",task.c_str(),etag.c_str()), Form("%s Parton #it{p}_{T}, #hat{#it{p}}_{T}=%.1f; #it{p}_{T} [GeV/#it{c}]", task.c_str(), pthat),100, 0., th1_ptmax));
      hgrams.push_back(new TH1D( Form("%s_parton_phi%s",task.c_str(),etag.c_str()), Form("%s Parton #phi, #hat{#it{p}}_{T}=%.1f; #phi [rad]", task.c_str(), pthat), 100, 0, 2*M_PI));
      hgrams.push_back(new TH1D( Form("%s_parton_eta%s",task.c_str(),etag.c_str()), Form("%s Parton #eta, #hat{#it{p}}_{T}=%.1f; #eta", task.c_str(), pthat), 100., -5., 5.));
      // add additional histograms as wanted; see class PythiaGun.h for options
    }
    break;

    case QA_TYPE::PREEQUILIBRIUM_DYNAMICS: {
      // preequilibrium dynamics just keeps std::vectors; get physical geometry from
      // the shaired pointer to InitialState for auto pre = std::dynamic_pointer_cast<PreequilibriumDynamics>(it->second.lock());
      auto pre = std::dynamic_pointer_cast<PreequilibriumDynamics>(it->second.lock());
      auto& istate = pre->ini;
      float xstep = istate->GetXStep();
      float xmax = istate->GetXMax();
      float ystep = istate->GetYStep();
      float ymax = istate->GetYMax();
      int nx = 2*xmax / xstep;
      int ny = 2*ymax / ystep;

      // See PreequiblibriumDynamics.h for available parameters
      // they are: e_, P_, utau_, ux_, uy_, ueta_, pi00_, 
      // pi01_, pi01_, pi02_, pi03_, pi11_, pi12_, pi13_, pi22_,
      // pi23_, pi33_, bulk_Pi,
      // get the 
      hgrams.push_back(new TH2D(Form("%s_e%s", task.c_str(), etag.c_str()), Form("%s Pre-Equilibrium Dynamics energy density [GeV/fm^{3}];x [fm]; y [fm]", task.c_str()), nx, -xmax, xmax, ny, -ymax, ymax));

      // FIXME: e_ is identical to Trento's entropy; why don't we use s?
    }
    break;

    case QA_TYPE::FLUID_DYNAMICS: {
      // Note: the information for the axes of fluid dynamics come from FluidEvolutionHistory,
      // which isn't available upon Init, therefore the histograms have to be initialized 
      // at the first Exec
      if (!has_first_exec) break;

      // NOTE: this information isn't availabe in the FluidEvolutionHistory upon Init,
      //       therefore the histograms have to be initialized upon the first Exec
      auto fluid = std::dynamic_pointer_cast<FluidDynamics>(it->second.lock());
      auto bInfo = fluid->get_bulk_info();

      int nX = bInfo.nx;
      double xMin = bInfo.XMin();
      double xMax = bInfo.XMax(); //same for y axis ...

      int nY = bInfo.ny;
      double yMin = bInfo.YMin();
      double yMax = bInfo.YMax(); //same for y axis ...

      // histograms can fill anything from FluidEvolutionHistory:
      // ENERGY_DENSITY, ENTROPY_DENSITY, TEMPERATURE,
      // PRESSURE, QGP_FRACTION, MU_B, MU_C, MU_S,
      // VX, VY, VZ, PI00, PI01, PI02,
      // PI03, PI11, PI12, PI13, PI22, PI23,
      // PI33, BULK_PI, INVALID
      // default here, just get the bulk entropy and temperature

      hgrams.push_back(new TH2D(Form("%s_entropy%s", task.c_str(), etag.c_str()), Form("%s initial entropy density [1/fm^{3}];x [fm]; y [fm]", task.c_str()), nX, -xMax, xMax, nY, -yMax, yMax));
      hgrams.push_back(new TH2D(Form("%s_T%s", task.c_str(), etag.c_str()), Form("%s initial temperature [GeV];x [fm]; y [fm]", task.c_str()), nX, -xMax, xMax, nY, -yMax, yMax));

      hgrams.push_back(new TH1D(Form("%s_CellFreezeOutTimes%s", task.c_str(), etag.c_str()), Form("%s fluid cells freezeout times; #tau [fm/c]", task.c_str()), 150, 0., 15.));
      hgrams.push_back(new TH1D(Form("%s_CellFreezeXlocs%s", task.c_str(), etag.c_str()), Form("%s fluid cells freezeout X locations; x [fm]", task.c_str()), 150, -xMax, xMax));
      hgrams.push_back(new TH1D(Form("%s_CellFreezeYlocs%s", task.c_str(), etag.c_str()), Form("%s fluid cells freezeout Y locations; y [fm]", task.c_str()), 150, -yMax, yMax));
      hgrams.push_back(new TH2D(Form("%s_FreezeOutEntropy%s", task.c_str(), etag.c_str()), Form("%s freezeout entropy [1/fm^{3}];x [fm]; y [fm]", task.c_str()), nX, -xMax, xMax, nY, -yMax, yMax));
    }
    break;

    case QA_TYPE::SOFT_PARTICLIZATION: {
      hgrams.push_back(new TH1D(Form("%s_Nhadrons%s", task.c_str(), etag.c_str()), Form("%s Number of Hadrons; N_{hadrons}", task.c_str()), th1_nhadrons, 0., th1_nhadrons));
      hgrams.push_back(new TH1D(Form("%s_hadron_pt%s", task.c_str(), etag.c_str()), Form("%s Hadron pT; #it{p}_{T} [GeV/#it{c}]", task.c_str()), 100, 0., th1_ptmax));
    }
    break;

    case QA_TYPE::JET_ENERGY_LOSS: {
      // get the initiating partons pTs of the jets
      hgrams.push_back(
          new TH1D(Form("%s_jet_init_parton_pt%s", task.c_str(), etag.c_str()),
                   Form("%s Jet Initiating Parton pT; #it{p}_{T} [GeV/#it{c}]",
                        task.c_str()),
                   100, 0., th1_ptmax));
      // get the jet pt
      hgrams.push_back(
          new TH1D(Form("%s_jet_pt%s", task.c_str(), etag.c_str()),
                   Form("%s R=0.7; #it{p}_{T} [GeV/#it{c}]", task.c_str()), 100,
                   0., th1_ptmax));
      // get the partons of the jet constituents
      hgrams.push_back(
          new TH1D(Form("%s_jet_const_pt%s", task.c_str(), etag.c_str()),
                   Form("%s R=0.7, Jet Constituent pT; #it{p}_{T} [GeV/#it{c}]",
                        task.c_str()),
                   100, 0., th1_ptmax));
      // get the jet constituent z distribution
      hgrams.push_back(new TH1D(
          Form("%s_jet_const_z%s", task.c_str(), etag.c_str()),
          Form("%s R=0.7, Jet Constituent z; z (p_{T}^{const.}/p_{T}^{jet})",
               task.c_str()),
          100, 0., 1.));
    } break;

    case QA_TYPE::HADRONIZATION: {
      // number of hadrons
      hgrams.push_back(
          new TH1D(Form("%s_Nhad%s", task.c_str(), etag.c_str()),
                   Form("%s Hadron; number of hadrons", task.c_str()),
                   th1_nhadrons, 0., 1. * th1_nhadrons));
      // hadron pT distribution
      hgrams.push_back(
          new TH1D(Form("%s_hadron_pt%s", task.c_str(), etag.c_str()),
                   Form("%s Hadron pT; #it{p}_{T} [GeV/#it{c}]", task.c_str()),
                   100, 0., th1_ptmax));

      // the rest of these mirror those of JET_ENERGY_LOSS...
      // get the jet pt
      hgrams.push_back(
          new TH1D(Form("%s_jet_pt%s", task.c_str(), etag.c_str()),
                   Form("%s R=0.7; #it{p}_{T} [GeV/#it{c}]", task.c_str()), 100,
                   0., th1_ptmax));
      // get the partons of the jet constituents
      hgrams.push_back(
          new TH1D(Form("%s_jet_const_pt%s", task.c_str(), etag.c_str()),
                   Form("%s R=0.7, Jet Constituent pT; #it{p}_{T} [GeV/#it{c}]",
                        task.c_str()),
                   100, 0., th1_ptmax));
      // get the jet constituent z distribution
      hgrams.push_back(new TH1D(
          Form("%s_jet_const_z%s", task.c_str(), etag.c_str()),
          Form("%s R=0.7, Jet Constituent z; z (p_{T}^{const.}/p_{T}^{jet})",
               task.c_str()),
          100, 0., 1.));
    } break;

    default:
    break;
  }
    return hgrams;
}

void JetScapeQA::FillHgrams(std::tuple<string, QA_TYPE, vector<TH1*>>& qa) {
  if (std::get<2>(qa).size()==0) {
    // this occurs when they are not generated in Init in the qa_tasks because
    // the FluidEvolutionHistory's geometry is not populated at time of Init
    std::get<2>(qa) = MakeHgrams(std::get<0>(qa), std::get<1>(qa));
  }

  auto it = taskMap.find(std::get<0>(qa));
  auto& qa_type = std::get<1>(qa);
  auto& hgrams = std::get<2>(qa);

  switch (qa_type) {
    case QA_TYPE::INITIAL_STATE: {
      auto iS = std::dynamic_pointer_cast<InitialState>(it->second.lock());

      // fill in centrality
      hgrams[1]->Fill(iS->GetEventCentrality());

      // fill in entropy
      auto entropy = iS->GetEntropyDensityDistribution();
      for (size_t index = 0; index < entropy.size(); ++index) {
        auto loc = iS->CoordFromIdx(index);
        static_cast<TH2D*>(hgrams[0])->Fill(std::get<0>(loc), std::get<1>(loc), entropy[index]);
      }
    }
    break;

    case QA_TYPE::HARD_PROCESS: {
      auto hp = std::dynamic_pointer_cast<HardProcess>(it->second.lock());
      auto inPartons = hp->GetPartonList();

      for (const auto& parton : inPartons) {
        hgrams[0]->Fill(parton->pt());
        hgrams[1]->Fill(parton->eta());
        hgrams[2]->Fill(parton->phi());
      }
    }
    break;

    case QA_TYPE::PREEQUILIBRIUM_DYNAMICS: {
      auto pre = std::dynamic_pointer_cast<PreequilibriumDynamics>(it->second.lock());
      auto istate = pre->ini;
      auto e = pre->e_;
      for (size_t index = 0; index < e.size(); ++index) {
        auto loc = istate->CoordFromIdx(index);
        static_cast<TH2D*>(hgrams[0])->Fill(std::get<0>(loc), std::get<1>(loc), e[index]);
      }
    }
    break;

    case QA_TYPE::FLUID_DYNAMICS: {
      auto fluid = std::dynamic_pointer_cast<FluidDynamics>(it->second.lock());
      auto bInfo = fluid->get_bulk_info();
      double tau0 = bInfo.Tau0();
      int nx = bInfo.nx;
      int ny = bInfo.ny;
      int id_eta = 0;

      // get the data for the initial distributions
      for (int ix = 0; ix < nx; ++ix) {
        auto x = bInfo.XCoord(ix);
        for (int iy=0; iy<ny;++iy) {
          auto y = bInfo.YCoord(iy);
          auto cell = bInfo.GetFluidCell(tau0, ix, iy, id_eta);
          static_cast<TH2D*>(hgrams[0])->Fill(x, y, cell.entropy_density);
          static_cast<TH2D*>(hgrams[1])->Fill(x, y, cell.temperature);
        }
      }

      // fill in freeze out information
      // this appears to be uniformly zero... :()

      std::vector<SurfaceCellInfo> m_surfaceCellVector;
      fluid->getSurfaceCellVector(m_surfaceCellVector);
      for (auto& cell : m_surfaceCellVector) {
        hgrams[2]->Fill(cell.tau);
        hgrams[3]->Fill(cell.x);
        hgrams[4]->Fill(cell.y);
        static_cast<TH2D*>(hgrams[5])->Fill(cell.x, cell.y, cell.entropy_density);
      }
    }
    break;

    case QA_TYPE::JET_ENERGY_LOSS: {
      // this task can be called mutiple times in the QA list -- only run once
      if (has_run_JEL) {
        has_run_JEL = true;
        break;
      }
        
      auto it = taskMap.equal_range("JetEnergyLoss");
      // get the initiating partons pTs of the jets
      for (auto itr = it.first; itr != it.second; ++itr) {
        auto init_Parton =
            std::dynamic_pointer_cast<JetEnergyLoss>(itr->second.lock())
                ->GetShowerInitiatingParton();
        hgrams[0]->Fill(init_Parton->pt());
      }

      // get the jet constituents and z
      vector<fjcore::PseudoJet> vfinals;
      // can do shower by shower QA here too, also same possible in the EbyE case ...
      for (auto itr = it.first; itr != it.second; ++itr)
      {
          auto mSfinal =  std::dynamic_pointer_cast<JetEnergyLoss>(itr->second.lock())->GetShower()->GetFinalPartonsForFastJet();
          vfinals.insert(vfinals.end(),mSfinal.begin(), mSfinal.end());     
      }
      fjcore::JetDefinition jet_def( fjcore::antikt_algorithm, 0.7); //hardcoded, make readable from XML maybe more R's ...
      fjcore::ClusterSequence hcs(vfinals, jet_def);
      vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets(min_jetpt));

      for (int k = 0; k < hjets.size(); k++) {
        auto jet_pt = hjets[k].pt();
        hgrams[1]->Fill(jet_pt);
        for (auto &c : hjets[k].constituents()) {
          hgrams[2]->Fill(c.pt());
          hgrams[3]->Fill(c.pt() / jet_pt); // z = pT_const
        }
      }
    }
    break;

    case QA_TYPE::HADRONIZATION: {
      auto hadro = std::dynamic_pointer_cast<Hadronization>(it->second.lock());

      hgrams[0]->Fill(hadro->GetHadrons().size()); // number of hadrons

      for (const auto& hadron : hadro->GetHadrons()) {
        hgrams[1]->Fill(hadron->pt());
      }

      vector<fjcore::PseudoJet> forFJ;
      for (auto &h : hadro->GetHadrons()) {
          forFJ.push_back(h->GetPseudoJet());
      }

      fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, 0.7); //hardcoded, make readable from XML maybe more R's ...
      fjcore::ClusterSequence hcs(forFJ, jet_def);
      vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets(min_jetpt));

      for (int k = 0; k < hjets.size(); k++) {
        auto jet_pt = hjets[k].pt();
        hgrams[2]->Fill(jet_pt);
        for (auto &c : hjets[k].constituents()) {
          hgrams[3]->Fill(c.pt());
          hgrams[4]->Fill(c.pt() / jet_pt); // z = pT_const
        }
      }
    }
    break;

    case QA_TYPE::SOFT_PARTICLIZATION: {
      auto soft = std::dynamic_pointer_cast<SoftParticlization>(it->second.lock());
      hgrams[0]->Fill(soft->Hadron_list_.size()); // number of hadrons
      for (auto& list : soft->Hadron_list_) {
        for (auto& hadron : list) {
          hgrams[1]->Fill(hadron->pt());
        }
      }
    }
    break;

    default:
        break;
  }
}

void JetScapeQA::Finish() {
    JSINFO << "Finish JetScapeQA : " << GetId() << " ...";

    //normalize per-event-histograms
    if (normalize_hgrams_per_event) {
      JSINFO << "Normalizing selected histograms per event ...";
      int nEvents = JetScapeModuleBase::GetCurrentEvent();
      for (auto& qa : qa_tasks) {
        for (auto& h : std::get<2>(qa)) {
          h->Scale(1.0 / (double) nEvents);
        }
      }
    }

    if (fOutputFile) {
      JSINFO<<"DEBUG OUTPUT FILE " << fOutputFile->GetName() << " is open, writing histograms ...";
        fOutputFile->Write();
        fOutputFile->Close();
        
        delete fOutputFile;
        fOutputFile = nullptr;
    }

    
    JSINFO << "JetScapeQA finished.";
    JSINFO << "JetScapeQA output file: " << outputFileName;

    PrintPDF();
}

void JetScapeQA::PrintPDF()
{
    JSINFO << "JetScapeQA::PrintPDF() to be implemented ...";
}


void JetScapeQA::UpdateTaskMap()
{
  VERBOSE(2) << "JetScapeQA::UpdateTaskMap()";

  //JP: Think about smarter/more efficient way rather than clear map and iterate through all tasks again ...
  taskMap.clear();
  auto mt = JetScapeSignalManager::Instance()->GetMainTaskPointer().lock();

  //Quick and dirty to see all tasks ... make recursive if needed
  if (mt) {
    for (auto it : mt->GetTaskList())
    {

      //JSINFO << t->GetId();
      taskMap.emplace(it->GetId(), it);

      for (auto it2 : it->GetTaskList())
      {
        //JSINFO  << it2->GetId() ;
        taskMap.emplace(it2->GetId(), it2);
      }
    }
  }
  PrintTaskMap();
}

vector<std::tuple<int,string,string>> JetScapeQA::GetTaskInfo(bool sorted) {
  // make a list of pairs of <taskName, taskDescription>
    std::vector<std::tuple<int,std::string, std::string>> taskInfo;
    for (auto &x : taskMap) {
      int taskNum = x.second.lock()->GetMyTaskNumber();
      std::stringstream ss;
      ss << " + " << std::left << std::setw(22) << x.first << ":"
         << x.second.lock().get()
         << "\t active = " << x.second.lock()->GetActive()
         << "\t Task number = " << taskNum;
      taskInfo.push_back({taskNum, x.first, ss.str()});
    }
    if (sorted) {
      std::sort(taskInfo.begin(), taskInfo.end());
    }
    return taskInfo;
}

void JetScapeQA::PrintTaskMap(bool sorted) {
  JSINFO << "JetScapeQA::PrintTaskMap()";
  auto taskInfo = GetTaskInfo(sorted);
  for (const auto &info : taskInfo) {
    JSINFO << std::get<2>(info);
  }
}

void JetScapeQA::PrintTasks()
{
  //Quick and dirty to see all tasks ... make recursive ...

  JSINFO << "JetScapeQA::PrintTasks()";

  auto mt = JetScapeSignalManager::Instance()->GetMainTaskPointer().lock();

  //Quick and dirty to see all tasks ... make recursive ...
  if (mt) {
    for (auto it : mt->GetTaskList()) {
      JSINFO << it->GetId();
      for (auto it2 : it->GetTaskList()) {
        JSINFO  << " " << it2->GetId() ;
        for (auto it3 : it2->GetTaskList())
          JSINFO  << "  " << it3->GetId() ;
      }
    }
  }
}

} // namespace Jetscape