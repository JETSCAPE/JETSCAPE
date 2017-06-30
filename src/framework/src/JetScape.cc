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

#include<iostream>

using namespace std;

namespace Jetscape {

JetScape::JetScape()
{
  n_events=1;
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
  
  INFO<<BOLDBLACK<<"Intialize JetScape ...";
  
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

  DEBUG<<"JetScape Debug from XML = "<< log_debug;
  DEBUG<<"JetScape Remark from XML = "<< log_remark;

   if ((int) m_vlevel>0)
     {
       JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
       DEBUG<<"JetScape Verbose Level from XML = "<<m_vlevel;
     }

   SetPointers();
   
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
      if (dynamic_pointer_cast<FluidDynamics>(it))
	JetScapeSignalManager::Instance()->SetHydroPointer(dynamic_pointer_cast<FluidDynamics>(it));
  
      if (dynamic_pointer_cast<JetEnergyLossManager>(it))
	JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(dynamic_pointer_cast<JetEnergyLossManager>(it));
      
      if (dynamic_pointer_cast<HardProcess>(it))
	JetScapeSignalManager::Instance()->SetHardProcessPointer(dynamic_pointer_cast<HardProcess>(it));
      
      if (dynamic_pointer_cast<JetScapeWriter>(it) && it->GetActive())
	JetScapeSignalManager::Instance()->SetWriterPointer(dynamic_pointer_cast<JetScapeWriter>(it)); 
    }
}

void JetScape::Exec()
{
  INFO<<BOLDBLACK<<"Run JetScape ...";
  INFO<<BOLDBLACK<<"Number of Events = "<<GetNumberOfEvents(); 
  
  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // --> JetScape is "Task Manager" of all modules ...

  // Simple way of passing the writer module pointer ...
  weak_ptr<JetScapeWriter> w;
  
  for (auto it : GetTaskList())
    if (dynamic_pointer_cast<JetScapeWriter>(it))
      if (it->GetActive())
	w=dynamic_pointer_cast<JetScapeWriter>(it);	           
  
  for (int i=0;i<GetNumberOfEvents();i++)
    {
      INFO<<BOLDBLACK<<"Run Event # = "<<i;
      DEBUG<<"Found "<<GetNumberOfTasks()<<" Modules Execute them ... ";
      
      JetScapeTask::ExecuteTasks();

      if (w.lock().get())
	JetScapeTask::WriteTasks(w);            
     
      JetScapeTask::ClearTasks();

      IncrementCurrentEvent();
    }
}

void JetScape::Finish()
{
  INFO<<BOLDBLACK<<"JetScape finished after "<<GetNumberOfEvents()<<" events!";
  DEBUG<<"More infos wrap up/saving to file/closing file ...";

  // same as in Init() and Exec() ...
  JetScapeTask::FinishTasks(); //dummy so far ...
}

} // end namespace Jetscape
