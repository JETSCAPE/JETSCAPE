// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScape.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "JetScapeBanner.h"
#include "InitialState.h"

#include<iostream>

using namespace std;

namespace Jetscape {

  /** Default constructor to create the main task of the JetScape framework. It sets the total number of events to 1.
   * By default, hydro events are used only once
   */
JetScape::JetScape() 
  : n_events(1)
  , reuse_hydro_ (false)
  , n_reuse_hydro_ (0)
{
  VERBOSE(8);
}

JetScape::~JetScape()
{
  VERBOSE(8);
  //JetScapeSignalManager::Instance()->Clear();
  //not needed, use weak_ptr in JetScapeSignalManager class (=not owning)
}

void JetScape::Show()
{
  show_jetscape_banner();
}

void JetScape::Init()
{
  Show();
  
  INFO<<BOLDRED<<"Intialize JetScape ...";
  
  JetScapeXML::Instance()->OpenXMLFile(GetXMLFileName());
  
  string log_debug = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "debug" )->GetText();
  string log_remark = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "remark" )->GetText();
  int m_vlevel=0;
  
  JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("vlevel")->QueryIntText(&m_vlevel);
  
  // Some general JetScape interface settings can also be read in from XML file ...
  // Order of task could also be verifed via XML or via sort of TaskList with order parameter
  if ((int) log_debug.find("off")>=0)
    JetScapeLogger::Instance()->SetDebug(false);
  
  if ((int) log_remark.find("on")>=0)
    JetScapeLogger::Instance()->SetRemark(true);

  JSDEBUG<<"JetScape Debug from XML = "<< log_debug;
  JSDEBUG<<"JetScape Remark from XML = "<< log_remark;

   if ((int) m_vlevel>0)
     {
       JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
       JSDEBUG<<"JetScape Verbose Level from XML = "<<m_vlevel;
     }

   SetPointers();
   
   // Set up helper. Mostly used for random numbers
   // Needs the XML reader singleton set up
   JSDEBUG<<"Seeding JetScapeTaskSupport from XML";
   JetScapeTaskSupport::ReadSeedFromXML( );

  // Has to be called explicitly since not really fully recursively (if ever needed)
   // So --> JetScape is "Task Manager" of all modules ...
   
   INFO<<"Found "<<GetNumberOfTasks()<<" Modules Initialize them ... ";
   JetScapeTask::InitTasks();
}

// kind of cluncky, maybe a better way ... ?
// Handle signal/slots in JetScape hence avoid passing pointers to sub tasks ...
void JetScape::SetPointers()
{
  
   // to get hydro pointer for signals, use signal?
  INFO<<"Set Hydro,JetEnergylossManager and IS Pointers for SignalManager to create Signal/Slots";
  
  for (auto it : GetTaskList())
    {
      if (dynamic_pointer_cast<InitialState>(it))
	JetScapeSignalManager::Instance()->SetInitialStatePointer(dynamic_pointer_cast<InitialState>(it));
 
      if (dynamic_pointer_cast<FluidDynamics>(it))
	JetScapeSignalManager::Instance()->SetHydroPointer(dynamic_pointer_cast<FluidDynamics>(it));
  
      if (dynamic_pointer_cast<JetEnergyLossManager>(it))
	JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(dynamic_pointer_cast<JetEnergyLossManager>(it));
      
      if (dynamic_pointer_cast<HardProcess>(it))
	JetScapeSignalManager::Instance()->SetHardProcessPointer(dynamic_pointer_cast<HardProcess>(it));
      
      if (dynamic_pointer_cast<JetScapeWriter>(it) && it->GetActive())
	JetScapeSignalManager::Instance()->SetWriterPointer(dynamic_pointer_cast<JetScapeWriter>(it)); 

      if (dynamic_pointer_cast<PartonPrinter>(it))
        JetScapeSignalManager::Instance()->SetPartonPrinterPointer(dynamic_pointer_cast<PartonPrinter>(it));
    }
}

void JetScape::Exec()
{
  INFO<<BOLDRED<<"Run JetScape ...";
  INFO<<BOLDRED<<"Number of Events = "<<GetNumberOfEvents();
  
  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // --> JetScape is "Task Manager" of all modules ...

  // Simple way of passing the writer module pointer ...
  weak_ptr<JetScapeWriter> w;
  //weak_ptr<PartonPrinter> p;  

  for (auto it : GetTaskList())
  {
    if (dynamic_pointer_cast<JetScapeWriter>(it))
    {  
      if (it->GetActive())
        w=dynamic_pointer_cast<JetScapeWriter>(it);	           
    }
  } 
 
  for (int i=0;i<GetNumberOfEvents();i++)
    {
      INFO<<BOLDRED<<"Run Event # = "<<i;
      JSDEBUG<<"Found "<<GetNumberOfTasks()<<" Modules Execute them ... ";
      
      JetScapeTask::ExecuteTasks();

      //JetScapeTask::GetPartons(p);

      if (w.lock().get())
	JetScapeTask::WriteTasks(w);            

      // For reusal, deactivate task after it has finished but before it gets cleaned up.
      if ( reuse_hydro_ ){
	if ( n_reuse_hydro_<=0 ){
	  WARN << " reuse_hydro is set, but n_reuse_hydro=" << n_reuse_hydro_;
	  throw std::runtime_error ("Incompatible reusal settings.");
	}
	for (auto it : GetTaskList()){
	  if ( ! dynamic_pointer_cast<FluidDynamics>(it)) continue;
	  if ( i%n_reuse_hydro_ == n_reuse_hydro_-1 ){
	    JSDEBUG << " i was " << i << " i%n_reuse_hydro_ = " << i%n_reuse_hydro_ << " --> ACTIVATING";
	    it->SetActive(true);
	  } else{
	    JSDEBUG << " i was " << i << " i%n_reuse_hydro_ = " << i%n_reuse_hydro_ << " --> DE-ACTIVATING";
	    it->SetActive(false);
	  }
	}
      }
	
      // Now clean up, only affects active taskjs
      JetScapeTask::ClearTasks();

      IncrementCurrentEvent();
    }
}

void JetScape::Finish()
{
  INFO<<BOLDBLACK<<"JetScape finished after "<<GetNumberOfEvents()<<" events!";
  JSDEBUG<<"More infos wrap up/saving to file/closing file ...";

  // same as in Init() and Exec() ...
  JetScapeTask::FinishTasks(); //dummy so far ...
}

} // end namespace Jetscape
