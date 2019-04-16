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
// JetScape Framework Run Macro
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

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
    
  // Logger settings (can be also set also via XML file, although note in that case they will apply only after they are initialized)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  JetScapeLogger::Instance()->SetVerboseLevel(0);
   
  Show();

  // Create main Jetscape task, and assign XML configuration files
  auto jetscape = make_shared<JetScape>();
  jetscape->SetXMLMasterFileName("../config/jetscape_master.xml");
  jetscape->SetXMLUserFileName("../config/jetscape_user.xml");

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
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------";
  INFO_NICE<<"| ... JetScape Framework ... |";
  INFO_NICE<<"------------------------------";
  INFO_NICE;
}
