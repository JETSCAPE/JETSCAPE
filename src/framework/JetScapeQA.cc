#include "JetScapeQA.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "JetScapeSignalManager.h"

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

    UpdateTaskMap();
    PrintTaskMap();
}

void JetScapeQA::Exec() {
    JSINFO << "Executing JetScapeQA : " << GetId() << " ...";

    UpdateTaskMap();
    PrintTasks();
    //PrintTaskMap();

    if (enableEbyEQA) {
        JSINFO << "Filling QA histograms per event ...";
        // Fill histograms or perform QA tasks here
        // Example: fOutputFile->cd(); // Change to the output file directory
        //         someHistogram->Fill(someValue);
    } else {
        JSINFO << "QA histograms for per event are disabled.";
    }
}

void JetScapeQA::Finish() {
    JSINFO << "Finish JetScapeQA : " << GetId() << " ...";
    
    if (fOutputFile) {

        fOutputFile->Write();
        fOutputFile->Close();
        
        delete fOutputFile;
        fOutputFile = nullptr;
    }
    
    JSINFO << "JetScapeQA finished.";
    JSINFO << "JetScapeQA output file: " << outputFileName;
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