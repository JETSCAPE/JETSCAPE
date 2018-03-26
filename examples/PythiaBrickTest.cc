// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// ------------------------------------------------------------
// JetScape Framework Brick Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterAscii.h"
//#include "JetScapeWriterAsciiGZ.h"
//#include "JetScapeWriterHepMC.h"

// User modules derived from jetscape framework clasess
// to be used to run Jetscape ...
#include "AdSCFT.h"
#include "Matter.h"
#include "ElossModulesTestMartini.h"
#include "brick_jetscape.h"
#include "Gubser_hydro_jetscape.h"
#include "PythiaGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "HadronizationModuleTest.h"
#include "ColorlessHad.h"

// Add initial state module for test
#include "TrentoInitial.h"

#include <chrono>
#include <thread>

using namespace std;

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);

  cout<<endl;
    
  // DEBUG=true by default and REMARK=false
  // can be also set also via XML file (at least partially)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  
  Show();

  // auto jetscape = make_shared<JetScape>("./jetscape_init_pythiagun.xml",10);
  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",200);
  jetscape->SetId("primary");
  // jetscape->set_reuse_hydro (true);
  // jetscape->set_n_reuse_hydro (10);

  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto trento = make_shared<TrentoInitial>();
  auto hydro = make_shared<Brick> ();
  //auto hydro = make_shared<GubserHydro> ();
  
  auto matter = make_shared<Matter> ();
  auto martini = make_shared<Martini> ();
  auto adscft = make_shared<AdSCFT> ();

  auto pythiaGun= make_shared<PythiaGun> ();

  auto printer = make_shared<PartonPrinter>();
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  auto hadroModule = make_shared<HadronizationModuleTest> ();
  auto colorless = make_shared<ColorlessHad> ();
  
  // only pure Ascii writer implemented and working with graph output ...
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  //auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");  
  //auto writer= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
  //writer->SetActive(false);

  //Remark: For now modules have to be added
  //in proper "workflow" order (can be defined via xml and sorted if necessary)  
  jetscape->Add(trento);
  jetscape->Add(pythiaGun);

  //Some modifications will be needed for reusing hydro events, so far
  //simple test hydros always executed "on the fly" ...
  jetscape->Add(hydro);

  jloss->Add(matter);
  // jloss->Add(martini);
  // jloss->Add(adscft);

  jlossmanager->Add(jloss);
  
  jetscape->Add(jlossmanager);

  jetscape->Add(printer);
  //hadro->Add(hadroModule);
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);


  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();
  
  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

  // "dummy" so far ...
  // Most thinkgs done in write and clear ...
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  cout<<endl;

// Some information is only known after the full run,
  // Therefore store information at the end of the file, in a footer
  writer->WriteComment ( "EVENT GENERATION INFORMATION" );
  Pythia8::Info& info = pythiaGun->info;
  std::ostringstream oss;
  oss.str(""); oss << "nTried    = " << info.nTried();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "nSelected = " << info.nSelected();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "nAccepted = " << info.nAccepted();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "sigmaGen  = " << info.sigmaGen();  
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "sigmaErr  = " << info.sigmaErr();
  writer->WriteComment ( oss.str() );

  oss.str(""); oss << "eCM  = " << info.eCM();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "pTHatMin  = " << pythiaGun->GetpTHatMin();
  writer->WriteComment ( oss.str() );
  oss.str(""); oss << "pTHatMax  = " << pythiaGun->GetpTHatMax();
  //writer->WriteComment ( oss.str() );
  //oss.str(""); oss << "JETSCAPE Random Seed  = " << JetScapeTaskSupport::Instance()->GetRandomSeed();
  writer->WriteComment ( oss.str() );
  writer->WriteComment ( "/EVENT GENERATION INFORMATION" );

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  // Print pythia statistics
  // pythiaGun->stat();

  // Demonstrate how to work with pythia statistics
  //Pythia8::Info& info = pythiaGun->info;
  cout << " nTried    = " << info.nTried() << endl;
  cout << " nSelected = " << info.nSelected()  << endl;
  cout << " nAccepted = " << info.nAccepted()  << endl;
  cout << " sigmaGen  = " <<   info.sigmaGen()  << endl;  
  cout << " sigmaErr  = " <<   info.sigmaErr()  << endl;
 
  
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"| Brick Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}
