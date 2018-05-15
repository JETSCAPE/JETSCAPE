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
// ------------------------------------------------------------
// JetScape Framework jet in hydro from file Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>
#include <thread>

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
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "HydroFromFile.h"
#include "PGun.h"
#include "PythiaGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"

#ifdef USE_HDF5
#include "InitialFromFile.h"
#endif
// using namespace std;
// Add initial state module for test
#include "TrentoInitial.h"

#include <chrono>
#include <thread>

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
  JetScapeLogger::Instance()->SetVerboseLevel(6);
   
  Show();

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml", 2000);
  // auto jetscape = make_shared<JetScape>("./jetscape_init_pythiagun.xml",5);
  jetscape->SetId("primary");
  jetscape->SetReuseHydro (true);
  jetscape->SetNReuseHydro (10001);
    
    auto pGun= make_shared<PGun> ();

  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto hydro = make_shared<HydroFromFile> ();
  //auto hydro = make_shared<GubserHydro> ();
  
  auto matter = make_shared<Matter> ();
  auto lbt = make_shared<LBT> ();
  auto martini = make_shared<Martini> ();
  auto adscft = make_shared<AdSCFT> ();
  //DBEUG: Remark:
  //does not matter unfortunately since not called recursively, done by JetEnergyLoss class ...
  //matter->SetActive(false);
  //martini->SetActive(false);
  // This works ... (check with above logic ...)
  //jloss->SetActive(false);

  //auto pythiaGun= make_shared<PythiaGun> ();

  auto printer = make_shared<PartonPrinter> ();

  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  auto hadroModule = make_shared<ColoredHadronization> ();
  auto colorless = make_shared<ColorlessHadronization> ();

  // only pure Ascii writer implemented and working with graph output ...
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  //auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");  
#ifdef USE_HEPMC
  //auto writer= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
#endif
  //writer->SetActive(false);

  //Remark: For now modules have to be added
  //in proper "workflow" order (can be defined via xml and sorted if necessary)
  
#ifdef USE_HDF5
  auto initial = make_shared<InitialFromFile>();
  jetscape->Add(initial);
#endif

//  jetscape->Add(pythiaGun);
    jetscape->Add(pGun);

   //Some modifications will be needed for reusing hydro events, so far
  //simple test hydros always executed "on the fly" ...
  jetscape->Add(hydro);

  // Matter with silly "toy shower (no physics)
  // and Martini dummy ...
  // Switching Q2 (or whatever variable used
  // hardcoded at 5 to be changed to xml)
  jloss->Add(matter);
  //jloss->Add(lbt);  // go to 3rd party and ./get_lbtTab before adding this module
  //jloss->Add(martini);
  //jloss->Add(adscft);
  
  jlossmanager->Add(jloss);
  
  jetscape->Add(jlossmanager);

  jetscape->Add(printer);

  hadro->Add(hadroModule);
  //hadro->Add(colorless);
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
/*  writer->WriteComment ( "EVENT GENERATION INFORMATION" );
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
*/
  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  // Print pythia statistics
  // pythiaGun->stat();

  // Demonstrate how to work with pythia statistics
  //Pythia8::Info& info = pythiaGun->info;
/*  cout << " nTried    = " << info.nTried() << endl;
  cout << " nSelected = " << info.nSelected()  << endl;
  cout << " nAccepted = " << info.nAccepted()  << endl;
  cout << " sigmaGen  = " <<   info.sigmaGen()  << endl;  
  cout << " sigmaErr  = " <<   info.sigmaErr()  << endl;
*/
  
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------------------------------";
  INFO_NICE<<"| Jet in hydro from file Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------------------------";
  INFO_NICE;
}
