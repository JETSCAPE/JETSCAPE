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
#include "JetScapeSignalManager.h"
#include "FluidDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "ElossModulesTest.h"
#include "sigslot.h"
 
using namespace std;
using namespace sigslot;

// -------------------------------------
// as shortcut if desired ...
// (can/should also be done in Classes if needed ...)

//#define Logger JetScapeLogger::Instance()
//#define INFO Logger->Info()<<__PRETTY_FUNCTION__
//#define VERBOSE(l) JetScapeLogger::Instance()->Verbose(l)<<__PRETTY_FUNCTION__
// Done and defined in JetScapeLogger.h ... 
// -------------------------------------

// Forward declaration

void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  cout<<endl;
  
  // DEBUG=true by default and REMARK=false
  // can be also set via XML file
  //JetScapeLogger::Instance()->SetDebug(false);
  //JetScapeLogger::Instance()->SetRemark(true);  
  JetScapeLogger::Instance()->SetVerboseLevel(9);
   
  Show();

  //JetScapeSignalManager::Instance();
  
  //Test verbose ...
  // Verbose Level 10 should be the highest ... (not check though in class) 
  //JetScapeLogger::Instance()->Verbose(4)<<"Test Verbose ...";
  //VERBOSE(4)<<" Test Verbose ...";
  
  //With definitions above/in logger class  as a shortcut
  //DEBUG<<" Test";
  
  //JetScapeXML::Instance()->SetXMLFileName("./jetscape_init.xml");
  //JetScapeXML::Instance()->OpenXMLFile();
  //JetScapeXML::Instance()->OpenXMLFile("./jetscape_init.xml");	
  
  //shared_ptr<JetScape> jetscape = make_shared<JetScape> ();
  // or shorter ...
  //auto jetscape = make_shared<JetScape> ();
  //jetscape->SetXMLFileName("./jetscape_init.xml");

  // JetScape Clas acts as TaskManager (should be generalized at some point)
  // Make it truly recursive (not yet really implemented that way)
  // Think harder and maybe in general/decide on shared vs. unique vs. raw pointer usage ...
  // Current: Just use make_shared always (propably not the most efficient solution ...)
  
  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",3);
  // if << overloaded for classes then for example (see commented out in JetScape.h)
  //cout<<*jetscape<<endl; 
  // try to automatically create in Add()!? and to handle things behind the scene ...
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();
  auto hydro = make_shared<FluidDynamics> ();
  
  auto matter = make_shared<Matter> ();
  auto martini = make_shared<Martini> ();

  jetscape->Add(hydro);

  jloss->Add(matter);
  jloss->Add(martini);

  //jetscape->Add(jloss);
  
  jlossmanager->Add(jloss);
  jetscape->Add(jlossmanager);
  
  // can add before or after jetscape (order does not matter)
  
  //INFO<<"Number of JetScape Tasks = "<<jetscape->GetNumberOfTasks();
  //INFO<<"Number of Eloss Tasks = "<<jloss->GetNumberOfTasks();
  
  jetscape->Init();

   // maybe smarter to deal with multiple eloss instances (think of di-jet) in seperate task Manager !??
  //auto jloss2=make_shared<JetEnergyLoss> (*jloss);
  //jetscape->Add(jloss2);

  REMARK<<"Module testing for now on Module Base class!";
  REMARK<<"jetscape->Init(); Calls all Init()'s of jetscape modules tasks, ";
  REMARK<<"since list/vector is sortable, can be used to ensure correct order or via xml file"; 

  // --------------------------
  /*
  //cout<<((JetEnergyLoss*) (jloss.get()))->GetQhat()<<endl;
  auto jloss2=make_shared<JetEnergyLoss> (*jloss); //default copy works with STL vector list but no deep copy!!? See below
  // to be implemented in class ...
  
  cout<<jloss->GetNumberOfTasks()<<" "<<jloss2->GetNumberOfTasks()<<endl;
  
  cout<<jloss.get()<<" "<<jloss2.get()<<endl;
  cout<<jloss->GetTaskAt(0).get()<<" "<<jloss2->GetTaskAt(0).get() <<endl;
  
  cout<<((Matter*) (jloss->GetTaskAt(0).get()))->GetQhat()<<endl; //check with shared ...
  cout<<((Matter*) (jloss2->GetTaskAt(0).get()))->GetQhat()<<endl;
  ((Matter*) (jloss2->GetTaskAt(0).get()))->SetQhat(12);
  cout<<((Matter*) (jloss->GetTaskAt(0).get()))->GetQhat()<<endl; //check with shared ...
  cout<<((Matter*) (jloss2->GetTaskAt(0).get()))->GetQhat()<<endl;
  */
  
  /*
  //cout<<((Martini*) (jloss->GetTaskAt(1).get()))->GetQhat()<<endl;
  //cout<<((Martini*) (jloss2->GetTaskAt(1).get()))->GetQhat()<<endl;
  */

  /*
  auto matter2 = make_shared<Matter> (*matter);
  matter2->SetQhat(15);
  cout<<matter2->GetQhat()<<endl;
  cout<<matter->GetQhat()<<endl;
  */
  // --------------------------
  
  //
  // Simple signal/slot creation
  //
  //hydro.get() --> old C++ style pointer needd for sigslot code ...
  //Maybe switch to boost signal2 for more modern implementation ... 
  //jloss->jetSignal.connect(hydro.get(), &FluidDynamics::UpdateEnergyDeposit);
  
  //Signal itself defined in JetEnergyLoss class
  //Signal not allowed virtual ... (think about inheritence and user interface ...!?)
  //matter->jetSignal.connect(hydro.get(), &FluidDynamics::UpdateEnergyDeposit);
  //matter2->jetSignal.connect(hydro.get(), &FluidDynamics::UpdateEnergyDeposit);
  //martini->jetSignal.connect(hydro.get(), &FluidDynamics::UpdateEnergyDeposit);

  //DEBUG<<"Connect Signal/Slot : jetSignal.connect(hydro.get(), &FluidDynamics::UpdateEnergyDeposit);";
  //REMARK<<"Can be handelded in own connection class ...";

  //matter2->jetSignal(1,2);
  
  // ------------------
  //some quick and dirty debug ...
  //cout<<jetscape.use_count()<<endl;
  //cout<<jloss.use_count()<<endl;
  
  //VERBOSE(9)<<" : Martini use count = "<<  martini.use_count() ;

  /*
  //martini = NULL;
  //martini.reset();
  //delete martini.get();
  JetScapeLogger::Instance()->Verbose(9)<<" Martini after nullptr ... "<<  martini.use_count() ; 

  cout<<jetscape.use_count()<<endl;
  cout<<jloss.use_count()<<endl;
  */
  // ------------------
  
  jetscape->Exec();
  
  //REMARK<<"Overload Exec for Hydro ...";

  // Others like finish and write (and in general data structure handling to be added/implmented accordingly)
  
  // If not handeled via task class can be called individually or put in jetscape main class ...
  //jloss->Init();
  //hydro->Init();
  
  //this_thread::sleep_for(std::chrono::milliseconds(1000000));
  
  INFO_NICE<<"Finished!";
  cout<<endl;
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-------------------------------";
  INFO_NICE<<"| Test JetScape Framework ... |";
  INFO_NICE<<"-------------------------------";
}
