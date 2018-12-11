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
// JetScape Framework hydro from file Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// Move it here to avoid conflicts:
// 1. Conflicts with MUSIC macros: hbarc, theta, limit, etc
// 2. Conflict of make_unique from smash and from JetScape
#include "SmashWrapper.h"

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "SoftHadronWriter.h"
//#ifdef USE_HEPMC
//#include "JetScapeWriterHepMC.h"
//#endif

// User modules derived from jetscape framework clasess
#include "AdSCFT.h"
#include "Matter.h"
#include "Martini.h"
#include "MusicWrapper.h"
// Make sure that nasty MUSIC macros are neutralized
#undef PI
#undef hbarc
#undef default_tol
#undef absol
#undef maxi
#undef mini
#undef sgn
#undef theta
#undef gmn
#undef limit

#include "iSpectraSamplerWrapper.h"
#include "TrentoInitial.h"
#include "PGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"

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
  //If you want to suppress it: use SetVerboseLevel(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(8);

  Show();

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",1);
  jetscape->SetReuseHydro (false);
  jetscape->SetNReuseHydro (0);

  // surface sampler
  auto iSS = make_shared<iSpectraSamplerWrapper> ();
  jetscape->Add(iSS);

  // afterburner
  auto smash = make_shared<SmashWrapper> ();
  jetscape->Add(smash);

  // Output
  auto writer= make_shared<SoftHadronWriterAscii> ("final_smash_hadrons.dat");
  // same as JetScapeWriterAscii but gzipped
  // auto writer= make_shared<SoftHadronWriterAsciiGZ> ("final_smash_hadrons.gz");
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

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-----------------------------------------------------";
  INFO_NICE<<"| Sampler_Afterburner :  JetScape Framework ... |";
  INFO_NICE<<"-----------------------------------------------------";
  INFO_NICE;
}
