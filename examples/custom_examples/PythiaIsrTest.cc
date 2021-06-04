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
// -------------------------------------------------
// XSCAPE Framework Clock Pythia Brick Test Program
// -------------------------------------------------

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
#include "BrickTest.h"
#include "GubserHydro.h"
#include "PythiaGun.h"
#include "InitialStateRadiationTest.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"
#include "CascadeTest.h"
#include "IsrManager.h"
#include "DummySplit.h"
#include "PartonShowerGeneratorDefault.h"
#include "IsrJet.h"
#include "IsrShowerPSG.h"

#include "MainClock.h"
#include "ModuleClock.h"

#include "QueryHistory.h"

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
  JetScapeLogger::Instance()->SetDebug(true);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);


  Show();

  // -------------
  //Test clock ...

  //auto mClock = make_shared<MainClock>();
  //mClock->SetTimeRefFrameId("SpaceTime");

  // clocks here are defaulted for testing, clocks can costumized via inhererting from the MainClock/ModuleClock base classes ...
  auto mClock = make_shared<MainClock>("SpaceTime",-1,3,0.1); // JP: make consistent with reading from XML in init phase ...
  auto mModuleClock = make_shared<ModuleClock>();
  mModuleClock->SetTimeRefFrameId("SpaceTime * 2");

  mClock->Info();
  mModuleClock->Info();

  auto jetscape = make_shared<JetScape>();
  jetscape->SetXMLMasterFileName("../config/jetscape_master.xml");
  jetscape->SetXMLUserFileName("../config/jetscape_user_test.xml");
  jetscape->SetId("primary");
  jetscape->AddMainClock(mClock);
  jetscape->ClockInfo();

  // Initial conditions and hydro
  //auto trento = make_shared<TrentoInitial>();
  auto trento = make_shared<InitialState>();
  auto pythiaGun= make_shared<PythiaGun> ();
  //auto isr = make_shared<InitialStateRadiationTest> ();
  auto hydro = make_shared<Brick> ();

  //auto hydroTest = make_shared<BrickTest> ();
  //hydroTest->SetMultiThread(true);
  //hydroTest->SetActive(false);

  jetscape->Add(trento);

  auto isrManager = make_shared<IsrManager>();
  //auto isrJloss = make_shared<JetEnergyLoss> (); //to be followed up ... make isr module ... !!!!
  auto isrJloss = make_shared<IsrJet>();
  auto oldPSG = make_shared<PartonShowerGeneratorDefault>(); //modify for ISR evolution ... to be discussed ...
  auto iDummy = make_shared<DummySplit> ();

  isrJloss->SetDeltaT(-0.1); isrJloss->SetStartT(0); isrJloss->SetMaxT(-3.); //will be moved to XML and proper Init() in IsrJet later ...

  //REMARK: Think a bit harder about directed graph creation and time direction !!!!! Graph inversion !???
  //make positve just for testing of iterating through a shower ...
  //isrJloss->SetDeltaT(0.1); isrJloss->SetStartT(0); isrJloss->SetMaxT(3);

  isrJloss->AddPartonShowerGenerator(oldPSG);
  isrJloss->Add(iDummy);

  isrManager->Add(isrJloss);

  pythiaGun->Add(isrManager);

  jetscape->Add(pythiaGun);
  //jetscape->Add(isr);
  jetscape->Add(hydro);
  //jetscape->Add(hydroTest);

  // Energy loss
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();

  //Set inactive task (per event) and with main clock attached do per time step for these modules ...
  //Needed to overwrite functions: CalculateTime() and ExecTime(), in these functions get
  //time, either main clock time or if module clock attached the tranformed time via: GetModuleCurrentTime();

  jlossmanager->SetActive(false);
  jloss->SetActive(false);

  //***************************************************************************
  //REMARK: Ordering of graph with negative times and iteration to be fixed!!!
  //        Invert graph ... certainly needs some more thinking ... !!!
  //        Otherwise seems to working fine, definitely for positive times!!!
  //***************************************************************************
  /*
  jlossmanager->SetUseIntialPartonShower(true);
  jloss->SetUseIntialPartonShower(true);

  auto isrPSG = make_shared<IsrShowerPSG>();
  jloss->AddPartonShowerGenerator(isrPSG);
  */

  //quick and dirty to check if module clock transformation is working conceptually ...
  //jloss->AddModuleClock(mModuleClock);

  //Matter is added but not executed, need to implement the per time step execution in JetEnergyLoss::DoShower()...
  auto matter = make_shared<Matter> ();
  auto dummy = make_shared<DummySplit> ();
  // auto lbt = make_shared<LBT> ();
  //auto martini = make_shared<Martini> ();
  //auto adscft = make_shared<AdSCFT> ();

  //Has to be set now if one wants to deal with negative
  //times in forward evolution; default is: 0 -- 100 ...
  jlossmanager->SetTimeRange(-20.0,20.0);
  jloss->SetTimeRange(-20.0,20.0);
  matter->SetTimeRange(-20.0,20.0);
  dummy->SetTimeRange(-20.0,20.0);

  // Note: if you use Matter, it MUST come first (to set virtuality)
  //jloss->Add(matter);
  jloss->Add(dummy);

  // jloss->Add(lbt);  // go to 3rd party and ./get_lbtTab before adding this module
  // jloss->Add(martini);
  //jloss->Add(adscft);
  jlossmanager->Add(jloss);
  jetscape->Add(jlossmanager);

  auto cascadeTest = make_shared<CascadeTest> ();
  cascadeTest->SetMultiThread(true);
  cascadeTest->SetActive(false);
  //jetscape->Add(cascadeTest);


  // JP: Leave out for now for testing clock(s) ... has to be updated accordingly ... (see JetEnergyLossManager as an example ...)
  // Hadronization
  /*
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  //auto hadroModule = make_shared<ColoredHadronization> ();
  //hadro->Add(hadroModule);
  auto colorless = make_shared<ColorlessHadronization> ();
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);
  */

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  writer->SetId("AsciiWriter"); //for task search test ...
  jetscape->Add(writer);

  /*
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
  */

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
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE<<"| Clock Brick Test XSCAPE Framework ...   |";
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE;
}
