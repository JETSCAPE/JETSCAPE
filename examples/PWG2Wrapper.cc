// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
// to be used to run Jetscape ...
// #include "AdSCFT.h"
#include "Matter.h"
// #include "Martini.h"
#include "Brick.h"
#include "PythiaGun.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"

// // Add initial state module for test
// #include "TrentoInitial.h"

#include <chrono>
#include <thread>
#include <sstream>

using namespace std;

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);
  
  string XMLname="./pwg2_init.xml";
  string outname="test_out.dat.gz";
  int Nevents = 10;  
  
  if ( argc >1 && string(argv[1]) == "-h" ) {
    cout << "Usage: PWG2Wrapper [xmlname] [outputname] [Nevents]" << endl;
    return -1;
  }

  switch ( argc ){
    break;
  case 4:
    Nevents=atoi( argv[3] );
    // Fallthrough
  case 3:
    outname = argv[2];
    // Fallthrough
  case 2:
    XMLname = argv[1];
    break;
  case 1:
    break;
  case 0:
    break;
  default:
    cout << "Usage: brickTest [xmlname] [outputname] [Nevents]" << endl;
    return -1;
  }
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

  auto jetscape = make_shared<JetScape>(XMLname,Nevents);
  jetscape->SetId("primary");
  // jetscape->set_reuse_hydro (true);
  // jetscape->set_n_reuse_hydro (10);

  // auto trento = make_shared<TrentoInitial>();
  auto init = make_shared<InitialState>();
  auto pythiaGun= make_shared<PythiaGun> ();
  auto hydro = make_shared<Brick> ();
  
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto matter = make_shared<Matter> ();

  // hadronization
  auto printer = make_shared<PartonPrinter>();
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  auto hadroModule = make_shared<ColoredHadronization> ();
  auto colorless = make_shared<ColorlessHadronization> ();
  
  // auto writer= make_shared<JetScapeWriterAscii> (outname);
  auto writer= make_shared<JetScapeWriterAsciiGZ> (outname);  
#ifdef USE_HEPMC
  // auto writer= make_shared<JetScapeWriterHepMC> (outname);
#endif

  //Remark: For now modules have to be added
  //in proper "workflow" order (can be defined via xml and sorted if necessary)  
  // jetscape->Add(trento);
  jetscape->Add(init);
  jetscape->Add(pythiaGun);
  jetscape->Add(hydro);

  // add module(s) to the eloss wrapper, than the eloss wrapper to the manager
  jloss->Add(matter);
  //jloss->Add(martini);
  //jloss->Add(adscft);
  jlossmanager->Add(jloss);  
  jetscape->Add(jlossmanager);

  // hadronization
  jetscape->Add(printer);
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);

  // Add the writer object
  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();
  
  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

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
  writer->WriteComment ( oss.str() );

  oss.str(""); oss << "JETSCAPE Random Seed  = " << JetScapeTaskSupport::Instance()->GetRandomSeed();
  writer->WriteComment ( oss.str() );

  writer->WriteComment ( "/EVENT GENERATION INFORMATION" );
  
  // Finalize
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

    
  // Print pythia statistics
  // pythiaGun->stat();

 
  
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
