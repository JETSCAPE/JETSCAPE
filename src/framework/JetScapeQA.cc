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

}