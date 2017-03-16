// ------------------------------------------------------------
// JetScape Framework Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <string>
#include <thread>

#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "ElossModulesTest.h"
#include "JetScapeWriterAsciiGZ.h"
#include "JetScapeWriterAscii.h"
//#include "JetScapeWriterHepMC.h"
//#include "HardProcess.h"

//#include "sigslot.h"

#include "brick_jetscape.h"
#include "Gubser_hydro_jetscape.h"
#include "PGun.h"
#include "JSPythia8.h"

#include "JetClass.hpp"

using namespace std;
using namespace sigslot;

// Forward declaration

void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  cout<<endl;

  // DEBUG=true by default and REMARK=false
  // can be also set via XML file
  JetScapeLogger::Instance()->SetDebug(true);
  JetScapeLogger::Instance()->SetRemark(true);  
  JetScapeLogger::Instance()->SetVerboseLevel(11);
   
  Show();

  //JetScapeSignalManager::Instance();

  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",3);
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto hydro = make_shared<Brick> ();
  //auto hydro = make_shared<GubserHydro> ();
  
  auto matter = make_shared<Matter> ();
  //auto matter2 = make_shared<Matter> ();
  auto martini = make_shared<Martini> ();

  auto pGun= make_shared<PGun> ();
  //auto py8=make_shared<JSPythia8> ("/Users/putschke/pythia8100/xmldoc",false);
  
  //auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  //auto writer= make_shared<JetScapeWriterHepMC> ("test_out.dat");
  //writer->SetActive(false);

  //jetscape->Add(py8);
  jetscape->Add(pGun);
  
  jetscape->Add(hydro);
  jloss->Add(matter);
  //jloss->Add(matter2);
  //jloss->Add(martini);
  
  jlossmanager->Add(jloss);
  jetscape->Add(jlossmanager);

  jetscape->Add(writer);
  
  jetscape->Init();
  
  jetscape->Exec();

  // "dummy" so far ...
  jetscape->Finish();
  
  INFO_NICE<<"Finished!";
  //cout<<endl;
   
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-------------------------------";
  INFO_NICE<<"| Test JetScape Framework ... |";
  INFO_NICE<<"-------------------------------";
}
