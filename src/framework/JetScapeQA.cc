#include "JetScapeQA.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "JetScapeSignalManager.h"
//#include "PythiaGun.h"

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<JetScapeQA> JetScapeQA::reg("JetScapeQA");

void JetScapeQA::Init() {
    JetScapeModuleBase::Init();
    JSINFO << "Initialize JetScapeQA : " << GetId() << " ...";
    
    enableEbyEQA = GetXMLElementInt({"JetScapeQA", "enableEbyEQA"});
    outputFileName = GetXMLElementText({"JetScapeQA", "outputFileName"});
    nEventsForQAHistograms = GetXMLElementInt({"JetScapeQA", "nEventsForQAHistograms"});

    JSINFO << "JetScapeQA: enableEbyEQA = " << enableEbyEQA;
    JSINFO << "JetScapeQA: outputFileName = " << outputFileName;
    JSINFO << "JetScapeQA: nEventsForQAHistograms = " << nEventsForQAHistograms;

    fOutputFile=new TFile(outputFileName.c_str(), "RECREATE");

    // read in from XML file too, binning etc ...
    hJetPartonPt = new TH1D("hJetPartonPt", "Jet Parton pT", 100, 0, 100);
    hJetHadronPt = new TH1D("hJetHadronPt", "Jet Hadron pT", 100, 0, 100);  

    UpdateTaskMap();
    PrintTaskMap();
}

void JetScapeQA::Exec() {
    VERBOSE(8) << "Executing JetScapeQA : " << GetId() << " ...";

    UpdateTaskMap();
    //PrintTasks();
    //PrintTaskMap();

    if (enableEbyEQA) {
        //JSINFO << "Filling QA histograms per event ...";
        // Fill histograms or perform QA tasks here
        // Example: fOutputFile->cd(); // Change to the output file directory
        //         someHistogram->Fill(someValue);
        if (GetCurrentEvent()<nEventsForQAHistograms)
            DoEbyEQA();
    }
    
    DoQA();
}

void JetScapeQA::DoEbyEQA() {
    JSINFO << "Performing event-by-event QA ...";
    // Implement the logic for event-by-event QA here
    // This could involve filling histograms, checking conditions, etc

    if (taskMap.count("PythiaGun")>0) HardProcessEbyEQA();
    if (taskMap.count("JetEnergyLoss")>0) JetPartonEbyEQA();

}   

void JetScapeQA::DoQA() {
    // JSINFO << "Performing general QA ...";
    // Implement the logic for general QA here
    // This could involve filling histograms, checking conditions, etc.
    if (taskMap.count("JetEnergyLoss")>0) JetPartonQA();
    if (taskMap.count("Hadronization")>0) JetHadronQA();
}

void JetScapeQA::HardProcessEbyEQA()
{
    string name = "InitalHardOutPt_event_" + std::to_string(GetCurrentEvent());
    TH1D *h = new TH1D(name.c_str(), "Hard Process EbyEQA Out Pt", 100, 0, 100);

    auto it = taskMap.find("PythiaGun");
    auto pyGun = std::dynamic_pointer_cast<HardProcess>(it->second.lock());

    auto inPartons = pyGun->GetPartonList();

    for (const auto& parton : inPartons) {
        h->Fill(parton->pt());
    }

    WriteEbyEQA(name, h);

    delete h;
}

void JetScapeQA::JetPartonEbyEQA()
{
    string name = "InitalHardInPt_event_" + std::to_string(GetCurrentEvent());
    TH1D *h = new TH1D(name.c_str(), "Jet Parton EbyEQA Hard In Pt", 100, 0, 100);

    int num = taskMap.count("JetEnergyLoss");
    auto it = taskMap.equal_range("JetEnergyLoss");

    for (auto itr = it.first; itr != it.second; ++itr)
    {
        auto inParton = std::dynamic_pointer_cast<JetEnergyLoss>(itr->second.lock())->GetShowerInitiatingParton();
        h->Fill(inParton->pt());
    }

    WriteEbyEQA(name, h);

    delete h;
}

void JetScapeQA::JetPartonQA()
{
    int num = taskMap.count("JetEnergyLoss");
    //DEBUG:
    //cout<<"--> "<<num<<endl;
    auto it = taskMap.equal_range("JetEnergyLoss");

    vector<fjcore::PseudoJet> vfinals;

    // can do shower by showe QA here too, also same possible in the EbyE case ...
    for (auto itr = it.first; itr != it.second; ++itr)
    {
        auto mSfinal =  std::dynamic_pointer_cast<JetEnergyLoss>(itr->second.lock())->GetShower()->GetFinalPartonsForFastJet();
        vfinals.insert(vfinals.end(),mSfinal.begin(), mSfinal.end());     
    }

    fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, 0.7); //hardcoded, make readable from XML maybe more R's ...
    fjcore::ClusterSequence hcs(vfinals, jet_def);
    vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets(2));

    for (int k=0;k<hjets.size();k++) {
	    //cout<<"Anti-kT jet "<<k<<" : "<<hjets[k].pt()<<endl;
        hJetPartonPt->Fill(hjets[k].pt());
    }
}

void JetScapeQA::JetHadronQA()
{
    auto it = taskMap.find("Hadronization");
    auto hadro = std::dynamic_pointer_cast<Hadronization>(it->second.lock());

    int nHadrons = hadro->GetHadrons().size();
    //cout<< "JetScapeQA::JetHadronQA() - Number of hadrons: " << nHadrons << endl;

    vector<fjcore::PseudoJet> forFJ;

    for (auto &h : hadro->GetHadrons()) {
        forFJ.push_back(h->GetPseudoJet());
    }
    
    //JP: Maybe make functio since in PartonQA same jet finding ...
    fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, 0.7); //hardcoded, make readable from XML maybe more R's ...
    fjcore::ClusterSequence hcs(forFJ, jet_def);
    vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets(2));

    for (int k=0;k<hjets.size();k++) {
	    //cout<<"Anti-kT jet "<<k<<" : "<<hjets[k].pt()<<endl;
        hJetHadronPt->Fill(hjets[k].pt());
    }

}

void JetScapeQA::SoftParticlizatonQA()
{

}

void JetScapeQA::Finish() {
    JSINFO << "Finish JetScapeQA : " << GetId() << " ...";

    NormalizePerEvent();

    if (fOutputFile) {

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

void JetScapeQA::NormalizePerEvent()
{
    JSINFO << "Normalizing selected histograms per event ...";

    int nEvents = JetScapeModuleBase::GetCurrentEvent();
    //cout<<nEvents<<endl;

    if (nEvents > 0) {
        hJetPartonPt->Scale(1.0 / (double) nEvents);
        hJetHadronPt->Scale(1.0 / (double) nEvents);
    }
}

void JetScapeQA::UpdateTaskMap()
{

  VERBOSE(2) << "JetScapeQA::UpdateTaskMap()";

  //JP: Think about smarter/more efficient way rather than clear map and iterate through all tasks again ...
  taskMap.clear();
  auto mt = JetScapeSignalManager::Instance()->GetMainTaskPointer().lock();
  //cout<<mt<<endl;

  //Quick and dirty to see all tasks ... make recursive if needed
  if (mt) {

    for (auto it : mt->GetTaskList())
    {

      //JSINFO << it->GetId();
      taskMap.emplace(it->GetId(), it);

      for (auto it2 : it->GetTaskList())
      {
        //JSINFO  << " " << it2->GetId() ;
        taskMap.emplace(it2->GetId(), it2);
      }
    }
  }
}

void JetScapeQA::PrintTaskMap()
{
  JSINFO << "JetScapeQA::PrintTaskMap()";

  for (auto& x : taskMap) {
    JSINFO << " " << x.first << ":\t " << x.second.lock().get() << "\t active = " << x.second.lock()->GetActive() << "\t Task number = "<<x.second.lock()->GetMyTaskNumber();
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

}