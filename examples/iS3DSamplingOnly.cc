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

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "iS3DWrapper.h"

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
  auto iS3D = make_shared<iS3DWrapper> ();
  jetscape->Add(iS3D);

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  // same as JetScapeWriterAscii but gzipped
  // auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  // HEPMC3
#ifdef USE_HEPMC
  // auto writer= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
#endif
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
  INFO_NICE<<"-----------------------------------------------";
  INFO_NICE<<"| iS3D Test JetScape Framework ... |";
  INFO_NICE<<"-----------------------------------------------";
  INFO_NICE;
}
