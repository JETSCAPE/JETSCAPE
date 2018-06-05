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
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBTD.h"
#include "Martini.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "PGun.h"
#include "PartonPrinter.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"

#include "tinyxml2.h"
#include <chrono>
#include <thread>

using namespace std;

using namespace Jetscape;

// Forward declaration
void Show();
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}


// -------------------------------------

void run_brick_test(double E0, double T, int events=3000)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);
  
  cout<<endl;
    
  // DEBUG=true by default and REMARK=false
  // can be also set also via XML file (at least partially)
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(8);
   
  Show();

  //modify the init.xml file
  JetScapeXML::Instance()->OpenXMLFile("./jetscape_init.xml");
  tinyxml2::XMLElement *pgunxml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PGun" )->FirstChildElement("E0"); 
  pgunxml->SetText(E0);
  tinyxml2::XMLElement *brickxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" )->FirstChildElement("Brick" )->FirstChildElement("T");
  brickxml->SetText(T);
  JetScapeXML::Instance()->CloseXMLFile();

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml", events);
  jetscape->SetId("primary");
  jetscape->SetReuseHydro(false);

/*
  JetScapeXML::Instance()->OpenXMLFile("./jetscape_init.xml");
  tinyxml2::XMLElement *pgunxml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PGun" );
  double pT=0.;
  pgunxml->FirstChildElement("pT")->QueryDoubleText(&pT);

  tinyxml2::XMLElement *brickxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" )->FirstChildElement("Brick" );
  double T=0.;
  brickxml->FirstChildElement("T")->QueryDoubleText(&T);
*/
  
  // Initial conditions and hydro
  auto trento = make_shared<TrentoInitial>();
  auto pGun= make_shared<PGun> ();
  auto hydro = make_shared<Brick> ();
  jetscape->Add(trento);
  jetscape->Add(pGun);
  jetscape->Add(hydro);

  // Energy loss
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto lbt = make_shared<LBTD> ();
  jloss->Add(lbt);
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
  std::string file_name="("+to_string_with_precision(E0,3)+", "+to_string_with_precision(T,3)+")test_out.dat";
  auto writer= make_shared<JetScapeWriterAscii> (file_name);
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
  
  INFO_NICE << file_name << "Finished!";
  cout<<endl;

  // wait for 5s
  //std::this_thread::sleep_for(std::chrono::milliseconds(500000));

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
  //printf ("Real time: %f seconds.\n",(start-end));
}

int main(int argc, char** argv)
{
  vector<double> E0vec;
  for(int i=0; i<14;i++)
  {
    double E0 = pow(10,0.7+0.1*i);
    E0vec.push_back(E0);
  }
  for(int i=1; i<10;i++)
  {
    E0vec.push_back(i*10.);
  }
  double E0;
  double T;
  for(int i=0; i<E0vec.size();i++)
  {
    for(int j=0;j<9;j++)
    {
      E0 = E0vec[i];
      T = 0.2+0.05*j;
      run_brick_test(E0,T);
    }
  }

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
