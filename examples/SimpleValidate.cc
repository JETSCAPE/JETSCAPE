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


#include <iostream>
#include <time.h>

// Add includes here to test if it breaks anything
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

#include "Brick.h"
#include "PGun.h"
#include "ElossValidation.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"

#include <chrono>

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);
  
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  JetScapeLogger::Instance()->SetVerboseLevel(0);
   
  Show();

  // Main framework task
  auto jetscape = make_shared<JetScape>("../examples/simplevalidate.xml",1);
  jetscape->SetId("primary");

  // Empty initial state
  auto ini = make_shared<InitialState>();
  ini->SetId("InitialState");

  // mono-energetic particle gun, fixed parameters in XML file
  auto pGun= make_shared<PGun> ();

  // Simple brick, parameters in XML file
  auto hydro = make_shared<Brick> ();

  // Energy loss manager, parameters in XML file
  auto jlossmanager = make_shared<JetEnergyLossManager> ();

  // Energy loss wrapper
  auto jloss = make_shared<JetEnergyLoss> ();
  
  // Energy loss module, can also add multiple ones.
  // Parameters in XML file
  // auto eloss1 = make_shared<ValidationEloss> ();
  auto eloss1 = make_shared<ElossValidate> ();
  
  // Pure Ascii writer
  auto writer= make_shared<JetScapeWriterAscii> ("validate_out.dat");
 
  //Remark: For now modules have to be added in proper "workflow" order
  jetscape->Add(ini);
  jetscape->Add(pGun);
  jetscape->Add(hydro);

  // add module to the eloss wrapper, than the eloss wrapper to the manager
  jloss->Add(eloss1);
  jlossmanager->Add(jloss);  
  jetscape->Add(jlossmanager);

  // TODO: Should add hadronizer here

  // Add the writer
  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

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
  INFO_NICE<<"--------------------------------------";
  INFO_NICE<<"| Validation Test JetScape Framework |";
  INFO_NICE<<"--------------------------------------";
  INFO_NICE;
}
