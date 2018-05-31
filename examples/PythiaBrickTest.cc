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
// JetScape Framework Brick Test Program with Pythia IS
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

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
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "PythiaGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"

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

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",200);
  jetscape->SetId("primary");

  // Initial conditions and hydro
  auto trento = make_shared<TrentoInitial>();
  auto pythiaGun= make_shared<PythiaGun> ();
  auto hydro = make_shared<Brick> ();
  jetscape->Add(trento);
  jetscape->Add(pythiaGun);
  jetscape->Add(hydro);


  // Energy loss
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();

  auto matter = make_shared<Matter> ();
  // auto lbt = make_shared<LBT> ();
  // auto martini = make_shared<Martini> ();
  // auto adscft = make_shared<AdSCFT> ();

  // Note: if you use Matter, it MUST come first (to set virtuality)
  jloss->Add(matter);
  // jloss->Add(lbt);  // go to 3rd party and ./get_lbtTab before adding this module
  // jloss->Add(martini);
  // jloss->Add(adscft);  
  jlossmanager->Add(jloss);  
  jetscape->Add(jlossmanager);

  
  // Hadronization
  // This helper module currently needs to be added for hadronization.
  auto printer = make_shared<PartonPrinter> ();
  jetscape->Add(printer);
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  auto hadroModule = make_shared<ColoredHadronization> ();
  hadro->Add(hadroModule);
  // auto colorless = make_shared<ColorlessHadronization> ();
  // hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);

  
  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  jetscape->Add(writer);
#ifdef USE_GZIP
  // same as JetScapeWriterAscii but gzipped
  auto writergz= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  jetscape->Add(writergz);
#endif
  // HEPMC3
#ifdef USE_HEPMC
  auto hepmcwriter= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
  jetscape->Add(hepmcwriter);
#endif

  // Intialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified
  jetscape->Exec();

  // For the future, cleanup is mostly already done in write and clear
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  // Print pythia statistics
  // pythiaGun->stat();

  // // Demonstrate how to work with pythia statistics
  // //Pythia8::Info& info = pythiaGun->info;
  // cout << " nTried    = " << info.nTried() << endl;
  // cout << " nSelected = " << info.nSelected()  << endl;
  // cout << " nAccepted = " << info.nAccepted()  << endl;
  // cout << " sigmaGen  = " <<   info.sigmaGen()  << endl;  
  // cout << " sigmaErr  = " <<   info.sigmaErr()  << endl;
   
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
